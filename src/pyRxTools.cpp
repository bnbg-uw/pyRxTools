/*
Copyright(C) 2024  University of Washington
This program is free software : you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.If not, see https ://www.gnu.org/licenses/.

Bryce Bartl - Geller
University of Washington Forest Resilience Lab
12 / 6 / 2024

cpprxgaming.cpp
*/

#include "pyRxTools.hpp"
#include "ReadProcessedFolder.hpp"
#include "allometry.hpp"
#include "vector.hpp"
#include "boost/stacktrace.hpp"
#include <chrono>
#include <ctime>

static std::vector<lapis::MultiPolygon> rxQ;

static bool lidarInit = FALSE;
static bool allomInit = FALSE;
static bool preProcessInit = FALSE;
static bool rxInit = FALSE;
std::unique_ptr<processedfolder::ProcessedFolder> lidarDataset;
static std::vector<RxGamingRxUnit> rxs;
static rxtools::allometry::UnivariateLinearModel allometry;

static std::default_random_engine dre;
static rxtools::Treatment treater;

static double convFactor;


void setProjDataDirectory(const char* searchpath) {
    proj_context_set_search_paths(nullptr, 1, &searchpath);

    auto x = proj_context_get_database_path(nullptr);
    std::cout << x << "\n";
}

void setSeed(int n) {
    dre = std::default_random_engine(n);
    treater = rxtools::Treatment(dre);
}

bool initLidarDataset(char* path) {
    try {
        lidarDataset = processedfolder::readProcessedFolder(path);
    }
    catch (std::exception e) {
        std::cerr << "Failed to read the lidar dataset. Is it properly formatted?\n";
        printf(e.what());
        return(FALSE);
    }
    lidarInit = TRUE;
    return(TRUE);
}

void reprojectPolygon(char* wkt, char* crsWkt, char* outWkt) {
    if (!lidarInit) {
        std::cout << "a" << "\n";
        throw std::runtime_error("Lidar dataset has not been initialized");
    }
    
    OGRGeometry* poGeom;
    auto err = OGRGeometryFactory::createFromWkt(&wkt, nullptr, &poGeom);
    if (err != OGRERR_NONE)
        throw std::invalid_argument("Failed to create geom from wkt");
    auto g = lapis::MultiPolygon(*poGeom, lapis::CoordRef{ crsWkt });
    OGRGeometryFactory::destroyGeometry(poGeom);

    auto outcrs = lidarDataset->crs();
    g.projectInPlace(outcrs);
    strcpy(outWkt, g.gdalGeometry()->exportToWkt().c_str());
}

bool queueRx(char* wkt, char* crsWkt) {
    if (lidarInit) {
        OGRGeometry* poGeom;
        auto err = OGRGeometryFactory::createFromWkt(&wkt, nullptr, &poGeom);
        if (err != OGRERR_NONE)
            throw std::invalid_argument("Failed to create geom from wkt");
        auto poly = lapis::MultiPolygon(*poGeom, lapis::CoordRef{ crsWkt });
        OGRGeometryFactory::destroyGeometry(poGeom);

        try {
            if (!poly.crs().isConsistent(lidarDataset->crs())) {
                poly.projectInPlace(lidarDataset->crs());
            }
        }
        catch (std::exception e) {
            std::cout << e.what();
            std::cout << "proj poly: " << crsWkt << "\n";
            std::cout << "proj lidar" << lidarDataset->crs().getPrettyWKT();

            throw e;
        }

        auto mask = lapis::Raster<int>(processedfolder::stringOrThrow(lidarDataset->maskRaster()));

        if (poly.overlaps(mask)) {
            rxQ.push_back(poly);
            return(true);
        }
        return(false);
    }
    else {
        throw std::runtime_error("Lidar dataset has not been initialized");
    }
}

void doPreProcessing(int nThread) {
    if (lidarInit) {
        std::cout << "Preprocessing:";
        std::pair<lapis::coord_t, lapis::coord_t> expectedRes{};
        auto units = lidarDataset->units();
        if (lidarDataset->type() == processedfolder::RunType::fusion) {
            if(units->name() == "metre") {
                expectedRes.first = 0.75;
            }
            else {
                expectedRes.first = 2.4606;
            }
        }
        else {
            expectedRes.first = lidarDataset->csmAlignment()->xres();
        }
        expectedRes.second = 0;

        std::cout << "RxQ size: " << rxQ.size() << "\n";
        rxs.resize(rxQ.size());

        int sofar = -1;
        std::mutex mut{};
        std::vector<std::thread> threads{};
        
        auto mwThreadFunc = [&](int i) {
            try {
                rastersAndTaosThread(nThread, i, mut, sofar, 2, 6, expectedRes);
            }
            catch (processedfolder::FileNotFoundException e) {
                std::cout << e.what() << '\n';
                std::cout << "Aborting\n";
                throw e;
            }
            catch (lapis::InvalidRasterFileException e) {
                std::cout << "Error reading file:\n";
                std::cout << e.what() << '\n';
                std::cout << "Aborting\b";
                throw e;
            }
            catch (std::exception e) {
                std::cout << e.what() << '\n';
                throw e;
            }
        };
        for (int i = 0; i < nThread; ++i) {
            threads.push_back(std::thread(mwThreadFunc, i));
        }
        for (int i = 0; i < nThread; ++i) {
            threads[i].join();
        }
    }
    else {
        throw std::runtime_error("Lidar dataset has not been initialized");
    }
    rxInit = TRUE;
}

void rastersAndTaosThread(const int nThread, const int thisThread, std::mutex& mut, int& sofar, const double canopycutoff,
    const double coregapdist, std::pair<lapis::coord_t, lapis::coord_t> expectedRes) {
    rxtools::TaoGettersPt getters(
        lapis::lico::alwaysAdd<lapis::VectorDataset<lapis::Point>>,
        lidarDataset->coordGetter(),
        [](const lapis::VectorDataset<lapis::Point>::ConstFeatureType& f) { return f.getRealField("height"); },
        [coregapdist](const lapis::VectorDataset<lapis::Point>::ConstFeatureType& f) { return coregapdist; },
        [](const lapis::VectorDataset<lapis::Point>::ConstFeatureType& f) { return 0.0; },
        [](const lapis::VectorDataset<lapis::Point>::ConstFeatureType& f) { return 0.0; }
    );
    //auto dist = 3 * lidarDataset.getConvFactor();
    //auto fixedRadius = [dist](lico::Tao& t) { t.crown = dist; };
    //areaproc::crownFunc cf = areaproc::crownFunc(fixedRadius);
    auto l = lapis::VectorDataset<lapis::MultiPolygon>(processedfolder::stringOrThrow(lidarDataset->tileLayoutVector()));
    auto mask = lapis::Raster<lapis::cell_t>(processedfolder::stringOrThrow(lidarDataset->maskRaster()));
    lapis::coord_t chmres = lapis::Raster<int>(processedfolder::stringOrThrow(lidarDataset->watershedSegmentRaster(0))).xres() * lidarDataset->getConvFactor();

    while (true) {
        mut.lock();
        ++sofar;
        if (sofar >= rxQ.size()) {
            mut.unlock();
            break;
        }
        int i = sofar;
        mut.unlock();
        auto poly = rxQ.at(i);
        std::cout << "Rx: " + std::to_string(sofar + 1) + "/" + std::to_string(rxQ.size()) + ". Starting on thread " + std::to_string(thisThread) + "\n";

        std::vector<lapis::Raster<int>> mhm;
        std::vector<lapis::Raster<double>> chm;
        std::vector<lapis::Raster<int>> basinMap;
        std::unordered_set<std::string> usedTiles;
        rxtools::TaoListPt taos;
        try {
            if (poly.overlaps(mask)) {
                for (int j = 0; j < l.nFeature(); j++) {
                    if (poly.overlaps(lidarDataset->extentByTile(j).value())) {
                        mut.lock();
                        auto thisMhm = lapis::Raster<int>(processedfolder::stringOrThrow(lidarDataset->maxHeightRaster(j)));
                        std::cout << "1";

                        mut.unlock();
                        std::cout << thisMhm.xres() << " " << thisMhm.yres() << "\n";
                        std::cout << expectedRes.first << " " << expectedRes.second << "\n";
                        //thisMhm.repairResolution(expectedRes.first, expectedRes.first, expectedRes.second, expectedRes.second);
                        std::cout << "2";

                        auto thisBasinMap = lapis::Raster<int>(processedfolder::stringOrThrow(lidarDataset->watershedSegmentRaster(j))) + i * 10000;
                        //thisBasinMap.repairResolution(expectedRes.first, expectedRes.first, expectedRes.second, expectedRes.second);
                        std::cout << "3";

                        thisMhm.defineCRS(mask.crs()); //is this line still needed?

                        mhm.push_back(thisMhm);
                        basinMap.push_back(thisBasinMap);

                        auto s = processedfolder::stringOrThrow(lidarDataset->csmRaster(j));
                        if (std::find(usedTiles.begin(), usedTiles.end(), s) == usedTiles.end()) {
                            auto thisChm = lapis::Raster<double>(s);
                            //thisChm.repairResolution(expectedRes.first, expectedRes.first, expectedRes.second, expectedRes.second);
                            chm.push_back(thisChm);
                            usedTiles.emplace(s);
                        }

                        rxtools::TaoListPt thisTaos{ lidarDataset->highPoints(j) };
                        for (int k = 0; k < thisTaos.size(); k++) {
                            if (poly.boundingBox().contains(thisTaos.x(k), thisTaos.y(k))) {
                                if (poly.containsPoint(lapis::Point(thisTaos.x(k), thisTaos.y(k), thisMhm.crs()))) {
                                    taos.addTAO(thisTaos[k]);
                                }
                            }
                        }
                    }
                } // for(int j = 0;...)

                auto mergeVectorInt = [](std::vector<lapis::Raster<int>>& v) {
                    std::vector<lapis::Raster<int>*> toMerge;
                    for (int i = 0; i < v.size(); i++) {
                        toMerge.push_back(&v[i]);
                    }
                    return lapis::mosaicInside(toMerge);
                };
                auto outMhm = mergeVectorInt(mhm);
                auto outBasin = mergeVectorInt(basinMap);

                auto mergeVectorDouble = [](std::vector<lapis::Raster<double>>& v) {
                    std::vector<lapis::Raster<double>*> toMerge;
                    for (int i = 0; i < v.size(); i++) {
                        toMerge.push_back(&v[i]);
                    }
                    return lapis::mosaicInside(toMerge);
                };
                auto outChm = mergeVectorDouble(chm);

                //This line may need to be adapted to the new code?
                //bigchm = spatial::resampleBilinear(bigchm, bigmhm);

                auto thisMask = mask;
                thisMask = lapis::cropRaster(thisMask, outMhm, lapis::SnapType::out);
                thisMask.maskByMultiPolygon(poly);
                thisMask = lapis::trimRaster(thisMask);

                outMhm.maskByMultiPolygon(poly);
                outChm.maskByMultiPolygon(poly);
                outBasin.maskByMultiPolygon(poly);

                outMhm = lapis::trimRaster(outMhm);
                outChm = lapis::trimRaster(outChm);
                outBasin = lapis::trimRaster(outBasin);

                outMhm = lapis::cropRaster(outMhm, outChm, lapis::SnapType::out);
                outMhm = lapis::extendRaster(outMhm, outChm, lapis::SnapType::in);
                outBasin = lapis::cropRaster(outBasin, outChm, lapis::SnapType::out);
                outBasin = lapis::extendRaster(outBasin, outChm, lapis::SnapType::in);

                auto rx = RxGamingRxUnit(thisMask, outMhm, outChm, outBasin, taos, convFactor);
                rxs.at(i) = rx;
            }
            else {
                continue;
            }
        } // try:
        catch (processedfolder::FileNotFoundException e) {
            mut.unlock();
            continue;
        }
    } //while (true)  #looping over rx units
}

int nTaos(int idx) {
    std::cout << "Rxs taos size: " << rxs[idx].taos.size() << "\n";
    return(rxs[idx].taos.size());
}

void getTaos(int idx, double* outData) {
    for (size_t i = 0; i < rxs[idx].taos.size(); ++i) {
        outData[5 * i] = rxs[idx].taos.x(i);
        outData[5 * i + 1] = rxs[idx].taos.y(i);
        outData[5 * i + 2] = rxs[idx].taos.area(i);
        outData[5 * i + 3] = rxs[idx].taos.height(i);
        outData[5 * i + 4] = rxs[idx].taos.radius(i);
    }
}

void setTaos(int idx, double* taoData, int size) {
    rxtools::TaoListPt taos;
    for (size_t i = 0; i < size / 5; ++i) {
        lico::Tao t = lico::Tao(
            taoData[5 * i], //x - 0
            taoData[5 * i + 1], //y - 1
            taoData[5 * i + 3], //height - 2
            taoData[5 * i + 4], //crown - 3
            taoData[5 * i + 2]); //area - 4
        taos.addTAO(t);
    }
    rxs[idx].taos = taos;
}

void setMhm(int idx, int* data, int nrow, int ncol, double xres, double yres, double xmin, double ymin, int naValue, char* projWKT) {
    auto a = lapis::Alignment(xmin, ymin, nrow, ncol, xres, yres, projWKT);
    rxs[idx].mhm = lapis::Raster<int>(a);
    for (lapis::cell_t c = 0; c < nrow*ncol; ++c) {
        rxs[idx].mhm[c].value() = data[c];
        if (data[c] == naValue || (naValue < -2000000 && data[c] < -2000000)) {
            rxs[idx].mhm[c].has_value() = false;
        }
        else {
            rxs[idx].mhm[c].has_value() = true;
        }
    }
}

void getMhmMeta(int idx, int* nrow, int* ncol, double* xres, double* yres, double* xmin, double* ymin, char* projWKT) {
    *nrow = rxs[idx].mhm.nrow();
    *ncol = rxs[idx].mhm.ncol();
    *xres = rxs[idx].mhm.xres();
    *yres = rxs[idx].mhm.yres();
    *xmin = rxs[idx].mhm.xmin();
    *ymin = rxs[idx].mhm.ymin();
    std::string wkt = rxs[idx].mhm.crs().getCompleteWKT();
    size_t bufferSize = wkt.length() + 1; //this is the value that wants to become the size of the buffer (the +1 is for the null terminator)
    strcpy_s(projWKT, bufferSize, wkt.c_str()); //using strcpy_s instead of strcpy makes it easy to not make mistakes, because it forces you to consider the buffer size
}

void getMhm(int idx, int* outData, int naValue) {
    for (lapis::cell_t i = 0; i < rxs[idx].mhm.ncell(); ++i) {
        if (rxs[idx].mhm[i].has_value()) {
            outData[i] = rxs[idx].mhm[i].value();
        }
        else {
            outData[i] = naValue;
        }
    }
};

void getChm(int idx, double* outData, double naValue) {
    for (lapis::cell_t i = 0; i < rxs[idx].chm.ncell(); ++i) {
        if (rxs[idx].chm[i].has_value()) {
            outData[i] = rxs[idx].chm[i].value();
        }
        else {
            outData[i] = naValue;
        }
    }
};


void setChm(int idx, double* data, int nrow, int ncol, double xres, double yres, double xmin, double ymin, double naValue, char* projWKT) {
    auto a = lapis::Alignment(xmin, ymin, nrow, ncol, xres, yres, projWKT);
    rxs[idx].chm = lapis::Raster<lapis::coord_t>(a);
    for (lapis::cell_t c = 0; c < nrow * ncol; ++c) {
        rxs[idx].chm[c].value() = data[c];
        if (data[c] == naValue || (naValue < -2000000 && data[c] < -2000000)) {
            rxs[idx].chm[c].has_value() = false;
        }
        else {
            rxs[idx].chm[c].has_value() = true;
        }
    }
}

void getChmMeta(int idx, int* nrow, int* ncol, double* xres, double* yres, double* xmin, double* ymin, char* projWKT) {
    *nrow = rxs[idx].chm.nrow();
    *ncol = rxs[idx].chm.ncol();
    *xres = rxs[idx].chm.xres();
    *yres = rxs[idx].chm.yres();
    *xmin = rxs[idx].chm.xmin();
    *ymin = rxs[idx].chm.ymin();
    std::string wkt = rxs[idx].chm.crs().getCompleteWKT();
    size_t bufferSize = wkt.length() + 1; //this is the value that wants to become the size of the buffer (the +1 is for the null terminator)
    strcpy_s(projWKT, bufferSize, wkt.c_str()); //using strcpy_s instead of strcpy makes it easy to not make mistakes, because it forces you to consider the buffer size
}


void setBasin(int idx, int* data, int nrow, int ncol, double xres, double yres, double xmin, double ymin, int naValue, char* projWKT) {
    auto a = lapis::Alignment(xmin, ymin, nrow, ncol, xres, yres, projWKT);
    rxs[idx].basinMap = lapis::Raster<int>(a);
    for (lapis::cell_t c = 0; c < nrow * ncol; ++c) {
        rxs[idx].basinMap[c].value() = data[c];
        if (data[c] == naValue || (naValue < -2000000 && data[c] < -2000000)) {
            rxs[idx].basinMap[c].has_value() = false;
        }
        else {
            rxs[idx].basinMap[c].has_value() = true;
        }
    }
}

void getBasinMeta(int idx, int* nrow, int* ncol, double* xres, double* yres, double* xmin, double* ymin, char* projWKT) {
    *nrow = rxs[idx].basinMap.nrow();
    *ncol = rxs[idx].basinMap.ncol();
    *xres = rxs[idx].basinMap.xres();
    *yres = rxs[idx].basinMap.yres();
    *xmin = rxs[idx].basinMap.xmin();
    *ymin = rxs[idx].basinMap.ymin();
    std::string wkt = rxs[idx].basinMap.crs().getCompleteWKT();
    size_t bufferSize = wkt.length() + 1; //this is the value that wants to become the size of the buffer (the +1 is for the null terminator)
    strcpy_s(projWKT, bufferSize, wkt.c_str()); //using strcpy_s instead of strcpy makes it easy to not make mistakes, because it forces you to consider the buffer size
}

void getBasin(int idx, int* outData, int naValue) {
    for (lapis::cell_t i = 0; i < rxs[idx].basinMap.ncell(); ++i) {
        if (rxs[idx].basinMap[i].has_value()) {
            outData[i] = rxs[idx].basinMap[i].value();
        }
        else {
            outData[i] = naValue;
        }
    }
};


void setMask(int idx, int* data, int nrow, int ncol, double xres, double yres, double xmin, double ymin, int naValue, char* projWKT) {
    auto a = lapis::Alignment(xmin, ymin, nrow, ncol, xres, yres, projWKT);
    rxs[idx].unitMask = lapis::Raster<lapis::cell_t>(a);
    rxs[idx].areaHa = 0;
    for (lapis::cell_t c = 0; c < nrow * ncol; ++c) {
        rxs[idx].unitMask[c].value() = data[c];
        if (data[c] == naValue || (naValue < -2000000 && data[c] < -2000000)) {
            rxs[idx].unitMask[c].has_value() = false;
        }
        else {
            rxs[idx].unitMask[c].has_value() = true;
            rxs[idx].areaHa += rxs[idx].unitMask.xres() * rxs[idx].unitMask.yres();
        }
    }
    rxs[idx].areaHa *= convFactor * convFactor;
    rxs[idx].areaHa /= 10000.0;
}

void getMaskMeta(int idx, int* nrow, int* ncol, double* xres, double* yres, double* xmin, double* ymin, char* projWKT) {
    *nrow = rxs[idx].unitMask.nrow();
    *ncol = rxs[idx].unitMask.ncol();
    *xres = rxs[idx].unitMask.xres();
    *yres = rxs[idx].unitMask.yres();
    *xmin = rxs[idx].unitMask.xmin();
    *ymin = rxs[idx].unitMask.ymin();
    std::string wkt = rxs[idx].unitMask.crs().getCompleteWKT();
    size_t bufferSize = wkt.length() + 1; //this is the value that wants to become the size of the buffer (the +1 is for the null terminator)
    strcpy_s(projWKT, bufferSize, wkt.c_str()); //using strcpy_s instead of strcpy makes it easy to not make mistakes, because it forces you to consider the buffer size
}

void getMask(int idx, int* outData, int naValue) {
    for (lapis::cell_t i = 0; i < rxs[idx].unitMask.ncell(); ++i) {
        if (rxs[idx].unitMask[i].has_value()) {
            outData[i] = rxs[idx].unitMask[i].value();
        }
        else {
            outData[i] = naValue;
        }
    }
};

void setAllometry(double intercept, double slope, int transform) {
    allometry.model.intercept = intercept;
    allometry.model.slope = slope;
    allometry.model.responseName = "DIA";
    // 0 = none
    // 1 = square
    // 2 = cube
    // 3 = log-log
    // 4 = suggest
    allometry.model.transform = (stats::Transform)transform;
    allometry.model.init = true;
    dbhFunc = [&](lico::adapt_type<spatial::unit_t> ht) {return allometry.getDbhFromHeightAuto(ht); };
    allomInit = TRUE;
}

void setAllometryFiaPath(char* path) {
    allometry.fiaPath = path;
}

int setAllometryWkt(char* wkt, char* crsWkt) {
    lapis::MultiPolygon poly;
    try {
        OGRGeometry* poGeom;
        auto err = OGRGeometryFactory::createFromWkt(&wkt, nullptr, &poGeom);
        if (err != OGRERR_NONE)
            throw std::invalid_argument("Failed to create geom from wkt");
        auto poly = lapis::MultiPolygon(*poGeom, lapis::CoordRef{ crsWkt });
        OGRGeometryFactory::destroyGeometry(poGeom);
    }
    catch (std::exception e) {
        std::cout << e.what() << "\n";
        std::cout << boost::stacktrace::stacktrace() << "\n";
        throw e;
    }

    auto allPlots = stats::FIALeastSquares::getPlotList(allometry.fiaPath);
    allometry.model = stats::FIALeastSquares(allPlots, allometry.fiaPath, "DIA");

    
    auto conv = lidarDataset->getConvFactor();
    std::cout << conv << "\n";
    lapis::Extent e(
        poly.boundingBox().xmin() - 10000 * conv,
        poly.boundingBox().xmax() + 10000 * conv,
        poly.boundingBox().ymin() - 10000 * conv,
        poly.boundingBox().ymax() + 10000 * conv,
        poly.crs());
    std::cout << poly.crs().getPrettyWKT() << "\n";
    auto n = allometry.model.limitByExtent(e);
    if (!n) {
        std::cout << "no fia plots in buffered extent of project area; aborting";
        throw std::runtime_error("no fia plots in buffered extent of project area");
    }

    allometry.model.initializeModel();
    dbhFunc = [&](lico::adapt_type<spatial::unit_t> ht) {return allometry.getDbhFromHeightAuto(ht); };

    allomInit = TRUE;
    return(n);
}

void getAllometry(double* intercept, double* slope, int* transform) {
    *intercept = allometry.model.intercept;
    *slope = allometry.model.slope;
    *transform = (int)allometry.model.transform;
}

//metric in, metric out
void getDbhFromHeight(double* height, double* dbh, int size) {
    for (int i = 0; i < size; ++i) {
        dbh[i] = allometry.getDbhFromHeightAuto(height[i]);
    }
}

void getCurrentStructure(int idx, double* ba, double* tph, double* mcs, double* cc) {
    *ba = rxs[idx].currentStructure.ba;
    *tph = rxs[idx].currentStructure.tph;
    *mcs = rxs[idx].currentStructure.mcs;
    *cc = rxs[idx].currentStructure.cc;
}

void calcCurrentStructure(int idx) {
    rxs[idx].currentStructure = rxtools::StructureSummary(rxs[idx].taos, align, rxs[idx].areaHa);
}

void setTargetStructure(int idx, double ba, double tph, double mcs, double cc) {
    rxs[idx].targetStructure = rxtools::StructureSummary(ba, tph, mcs, cc);
}

void getTargetStructure(int idx, double* ba, double* tph, double* mcs, double* cc) {
    *ba = rxs[idx].targetStructure.ba;
    *tph = rxs[idx].targetStructure.tph;
    *mcs = rxs[idx].targetStructure.mcs;
    *cc = rxs[idx].targetStructure.cc;
}

void getRawClumps(int idx, int* ids) {
    auto nd = lico::SparseDistMatrix();
    nd.nearDist(rxs[idx].taos.x(), rxs[idx].taos.y(), rxs[idx].taos.crown(), false);
    auto cl = nd.getClumps(rxs[idx].taos.size());
    for (size_t i = 0; i < cl.clumpID.size(); ++i) {
        ids[i] = cl.clumpID[i];
    }
}

void makeClumpMap(int idx, int* groupsizes, int* outData, int naValue) {
    try {
        auto b = rxs[idx].basinMap;

        std::unordered_map<int, int> taoIds;

        for (size_t i = 0; i < rxs[idx].taos.size(); ++i) {
            auto e = b.extract(rxs[idx].taos.x(i), rxs[idx].taos.y(i));
            if (e.has_value() && e.value() != 1) {
                taoIds.emplace(std::make_pair(e.value(), groupsizes[i]));
            }
        }

        for (lapis::cell_t j = 0; j < b.ncell(); ++j) {
            if (b[j].has_value()) {
                auto x = taoIds.find(b[j].value());
                if (x != taoIds.end()) {
                    b[j].value() = x->second;
                }
                else {
                    b[j].value() = 0;
                }
            }
        }

        for (lapis::cell_t i = 0; i < b.ncell(); ++i) {
            if (b[i].has_value()) {
                outData[i] = b[i].value();
            }
            else {
                outData[i] = naValue;
            }
        }
    }
    catch (std::exception e) {
        std::cout << e.what() << "\n";
        throw e;
    }
}

void getSimulatedStructures(int idx, double bbDbh, double* out) {
    try {
        auto rx = rxs[idx];

        std::vector<rxtools::StructureSummary> structures;

        auto dbh = dbhFunc(rx.taos.height());
        std::deque<size_t> notBBbase;
        lico::TaoList base;
        for (size_t i = 0; i < rx.taos.size(); i++) {
            if (dbh[i] >= bbDbh) {
                base.addTAO(rx.taos[i]);
            }
            else
            {
                notBBbase.push_back(i);
            }
        }
        auto testTaos = base;
        auto notBBidx = notBBbase;

        //double osi = getOsi(testTaos);
        double osi = 0;
        structures.push_back(rxtools::StructureSummary(testTaos, align, rx.areaHa));
        int step = (rx.taos.size() - testTaos.size()) / 10;

        //short to tall.
        for (int i = 1; i < 11; ++i) {
            auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            std::cout << i << " " << ctime(&t) << "\n";
            int limit = std::min<int>(rx.taos.size(), testTaos.size() + step);
            while (testTaos.size() < limit && notBBidx.size()) {
                auto idx = notBBidx.back();
                notBBidx.pop_back();
                testTaos.addTAO(rx.taos[idx]);
            }
            //osi = getOsi(testTaos);
            osi = 1;
            structures.push_back(rxtools::StructureSummary(testTaos, align, rx.areaHa));
        }

        //tall to short.
        testTaos = base;
        notBBidx = notBBbase;
        for (int i = 1; i < 11; ++i) {
            auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            int limit = std::min<int>(rx.taos.size(), testTaos.size() + step);
            while (testTaos.size() < limit && notBBidx.size()) {
                auto idx = notBBidx.front();
                notBBidx.pop_front();
                testTaos.addTAO(rx.taos[idx]);
            }
            //osi = getOsi(testTaos);
            osi = 2;
            structures.push_back(rxtools::StructureSummary(testTaos, align, rx.areaHa));
        }

        //random
        testTaos = base;
        notBBidx = notBBbase;
        std::shuffle(std::begin(notBBidx), std::end(notBBidx), dre);
        for (int i = 1; i < 11; ++i) {
            auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            int limit = std::min<int>(rx.taos.size(), testTaos.size() + step);
            while (testTaos.size() < limit && notBBidx.size()) {
                auto idx = notBBidx.back();
                notBBidx.pop_back();
                testTaos.addTAO(rx.taos[idx]);
            }
            //osi = getOsi(testTaos);
            osi = 3;
            structures.push_back(rxtools::StructureSummary(testTaos, align, rx.areaHa));
        }

        for (int i = 0; i < structures.size(); ++i) {
            for (int j = 0; j < 5; ++j) {
                out[5 * i + j] = structures[i][j];
                //if (j == 0) std::cout << out[4 * i + j] << "\n";
            }
        }
    }
    catch (std::exception e) {
        std::cout << e.what() << "\n";
        throw e;
    }
}

void doTreatment(int idx, double dbhMin, double dbhMax) {
    if (!rxInit) throw std::runtime_error("Rx not init");
    try {
        auto trt =  treater.doTreatmentGraph(rxs[idx], dbhMin, dbhMax, 10, false, "E:/dropbox/rxgaming paper/treatmentvisual/data/");
        rxs[idx].treatedTaos = std::get<0>(trt);
        rxs[idx].cutTaos = std::get<1>(trt);
        rxs[idx].result = std::get<2>(trt);
        rxs[idx].treatedStructure = rxs[idx].summarizeStructure(rxs[idx].treatedTaos, 0, dbhFunc);
    }
    catch (std::exception e) {
        std::cout << e.what();
        throw(e);
    }
}

void getTreatedChm(int idx, double* outData, double naValue) {
    auto rx = rxs[idx];
    std::unordered_set<int> basinIds;
    for (size_t i = 0; i < rx.treatedTaos.size(); ++i) {
        auto v = rx.basinMap.extract(rx.treatedTaos.x(i), rx.treatedTaos.y(i));
        if (v.has_value()) {
            basinIds.emplace(v.value());
        }
    }

    auto thisChm = rx.chm;
    for (lapis::cell_t i = 0; i < thisChm.ncell(); ++i) {
        if (thisChm[i].has_value()) {
            if (basinIds.find(rx.basinMap[i].value()) == basinIds.end()) {
                thisChm[i].value() = 0;
            }
        }
    }

    for (lapis::cell_t i = 0; i < thisChm.ncell(); ++i) {
        if (thisChm[i].has_value()) {
            outData[i] = thisChm[i].value();
        }
        else {
            outData[i] = naValue;
        }
    }
}

void getTreatedBasin(int idx, int* outData, int naValue) {
    auto rx = rxs[idx];
    std::unordered_set<int> basinIds;
    for (size_t i = 0; i < rx.treatedTaos.size(); ++i) {
        auto v = rx.basinMap.extract(rx.treatedTaos.x(i), rx.treatedTaos.y(i), lapis::ExtractMethod::near);
        if (v.has_value()) {
            basinIds.emplace(v.value());
        }
    }

    auto thisBasin = rx.basinMap;
    for (lapis::cell_t i = 0; i < thisBasin.ncell(); ++i) {
        if (thisBasin[i].has_value()) {
            if (basinIds.find(rx.basinMap[i].value()) == basinIds.end()) {
                thisBasin[i].value() = 1;
            }
        }
    }
    for (lapis::cell_t i = 0; i < thisBasin.ncell(); ++i) {
        if (thisBasin[i].has_value()) {
            outData[i] = thisBasin[i].value();
        }
        else {
            outData[i] = naValue;
        }
    }
}

int getNTreatedTaos(int idx) {
    return rxs[idx].treatedTaos.size();
}

void getTreatedTaos(int idx, double* outData) {
    for (size_t i = 0; i < rxs[idx].treatedTaos.size(); ++i) {
        outData[5 * i] = rxs[idx].treatedTaos.x(i);
        outData[5 * i + 1] = rxs[idx].treatedTaos.y(i);
        outData[5 * i + 2] = rxs[idx].treatedTaos.area(i);
        outData[5 * i + 3] = rxs[idx].treatedTaos.height(i);
        outData[5 * i + 4] = rxs[idx].treatedTaos.radius(i);
    }
}

int getNCutTaos(int idx) {
    return rxs[idx].cutTaos.size();
}

void getCutTaos(int idx, double* outData) {
    for (size_t i = 0; i < rxs[idx].cutTaos.size(); ++i) {
        outData[5 * i] = rxs[idx].cutTaos.x(i);
        outData[5 * i + 1] = rxs[idx].cutTaos.y(i);
        outData[5 * i + 2] = rxs[idx].cutTaos.area(i);
        outData[5 * i + 3] = rxs[idx].cutTaos.height(i);
        outData[5 * i + 4] = rxs[idx].cutTaos.radius(i);
    }
}

int getTreatmentResult(int idx) {
    return((int)rxs[idx].result);
}

void getTreatedStructure(int idx, double* ba, double* tph, double* mcs, double* cc) {
    *ba = rxs[idx].treatedStructure.ba;
    *tph = rxs[idx].treatedStructure.tph;
    *mcs = rxs[idx].treatedStructure.mcs;
    *cc = rxs[idx].treatedStructure.cc;
}

void getTreatedRawClumps(int idx, int* ids) {
    auto nd = lico::SparseDistMatrix();
    nd.nearDist(rxs[idx].treatedTaos.x(), rxs[idx].treatedTaos.y(), 6);
    auto cl = nd.getClumps(rxs[idx].taos.size());
    for (lico::index_t i = 0; i < cl.clumpID.size(); ++i) {
        ids[i] = cl.clumpID[i];
    }
}

void getTreatedClumpMap(int idx, int* inData, int inNaValue, int* groupsizes, int* outData, int naValue) {
    try {
        auto rx = rxs[idx];
        std::unordered_set<int> basinIds;
        for (size_t i = 0; i < rx.treatedTaos.size(); ++i) {
            auto v = rx.basinMap.extract(rx.treatedTaos.x(i), rx.treatedTaos.y(i));
            if (v.has_value()) {
                basinIds.emplace(v.value());
            }
        }
        auto b = rx.basinMap;
        for (lapis::cell_t c = 0; c < b.ncell(); ++c) {
            b[c].value() = inData[c];
            if (inData[c] == inNaValue || (inNaValue < -2000000 && inData[c] < -2000000)) {
                b[c].has_value() = false;
            }
            else {
                b[c].has_value() = true;
            }
        }

        std::unordered_map<int, int> taoIds;

        for (size_t i = 0; i < rx.treatedTaos.size(); ++i) {
            auto e = b.extract(rx.treatedTaos.x(i), rx.treatedTaos.y(i), lapis::ExtractMethod::near);
            if (e.has_value() && e.value() != 1) {
                taoIds.emplace(std::make_pair(e.value(), groupsizes[i]));
            }
        }

        for (lapis::cell_t j = 0; j < b.ncell(); ++j) {
            if (b[j].has_value()) {
                auto x = taoIds.find(b[j].value());
                if (x != taoIds.end()) {
                    b[j].value() = x->second;
                }
                else {
                    b[j].value() = 0;
                }
            }
        }

        for (lapis::cell_t i = 0; i < b.ncell(); ++i) {
            if (b[i].has_value()) {
                outData[i] = b[i].value();
            }
            else {
                outData[i] = naValue;
            }
        }

    }
    catch (std::exception e) {
        std::cout << e.what() << "\n";
        throw e;
    }
}

void setRxsSize(int i) {
    rxs.resize(i);
    rxInit = TRUE;
}

void setConvFactor(double cf) {
    convFactor = cf;
}

void getConvFactor(double* cf) {
    *cf = convFactor;
}