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

static std::vector<lapis::MultiPolygon> rxQ;

static bool lidarInit = FALSE;
static bool allomInit = FALSE;
static bool preProcessInit = FALSE;
static bool rxInit = FALSE;
std::unique_ptr<processedfolder::ProcessedFolder> lidarDataset;
static rxtools::TaoGettersPt getters;
static std::vector<RxGamingRxUnit> rxs;
static rxtools::allometry::UnivariateLinearModel dbhModel;

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
        getters = rxtools::TaoGettersPt(
            lapis::lico::alwaysAdd<lapis::VectorDataset<lapis::Point>>,
            lidarDataset->coordGetter(),
            lidarDataset->heightGetter(),
            lapis::lico::FixedRadius<lapis::VectorDataset<lapis::Point>>(3),
            lidarDataset->areaGetter(),
            [](const lapis::VectorDataset<lapis::Point>::ConstFeatureType& f)->double { throw std::runtime_error("dbhgetternotinit"); }
        );
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
        fprintf(stderr, "Lidar dataset has not been initialized");
        std::abort();
    }
    
    OGRGeometry* poGeom;
    auto err = OGRGeometryFactory::createFromWkt(&wkt, nullptr, &poGeom);
    if (err != OGRERR_NONE) {
        fprintf(stderr, "Failed to create geom from wkt");
        std::abort();
    }
    auto g = lapis::MultiPolygon(*poGeom, lapis::CoordRef{ crsWkt });
    OGRGeometryFactory::destroyGeometry(poGeom);

    auto outcrs = lidarDataset->crs();
    g.projectInPlace(outcrs);
    size_t bufferSize = g.gdalGeometry()->exportToWkt().length() + 1; //this is the value that wants to become the size of the buffer (the +1 is for the null terminator)
    strcpy_s(outWkt, bufferSize, g.gdalGeometry()->exportToWkt().c_str());
}

bool queueRx(char* wkt, char* crsWkt) {
    if (lidarInit) {
        OGRGeometry* poGeom;
        auto err = OGRGeometryFactory::createFromWkt(&wkt, nullptr, &poGeom);
        if (err != OGRERR_NONE) {
            fprintf(stderr, "Failed to create geom from wkt");
            std::abort();
        }
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

            std::abort();
        }

        auto mask = lapis::Raster<int>(processedfolder::stringOrThrow(lidarDataset->maskRaster()));

        if (mask.dataOverlapsMultiPolygon(poly)) {
            rxQ.push_back(poly);
            return(true);
        }
        return(false);
    }
    else {
        fprintf(stderr, "Lidar dataset has not been initialized");
        std::abort();
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
                std::abort();
            }
            catch (lapis::InvalidRasterFileException e) {
                std::cout << "Error reading file:\n";
                std::cout << e.what() << '\n';
                std::cout << "Aborting\b";
                std::abort();
            }
            catch (std::exception e) {
                std::cout << e.what() << '\n';
                std::abort();
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
        fprintf(stderr, "Lidar dataset has not been initialized");
        std::abort();
    }
    rxInit = TRUE;
}

void rastersAndTaosThread(const int nThread, const int thisThread, std::mutex& mut, int& sofar, const double canopycutoff,
    const double coregapdist, std::pair<lapis::coord_t, lapis::coord_t> expectedRes) {
    
    //auto dist = 3 * lidarDataset.getConvFactor();
    //auto fixedRadius = [dist](lico::Tao& t) { t.crown = dist; };
    //areaproc::crownFunc cf = areaproc::crownFunc(fixedRadius);
    auto l = lapis::VectorDataset<lapis::MultiPolygon>(processedfolder::stringOrThrow(lidarDataset->tileLayoutVector()));
    auto mask = lapis::Raster<lapis::cell_t>(processedfolder::stringOrThrow(lidarDataset->maskRaster()));

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
            if (mask.dataOverlapsMultiPolygon(poly)) {
                for (int j = 0; j < l.nFeature(); j++) {
                    if (poly.overlapsExtent(lidarDataset->extentByTile(j).value())) {
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

                        rxtools::TaoListPt thisTaos{ lidarDataset->highPoints(j)->string(), getters };
                        auto polyE = poly.boundingBox();
                        for (int k = 0; k < thisTaos.size(); k++) {
                            if (polyE.contains(thisTaos.x(k), thisTaos.y(k))) {
                                if (poly.containsPoint(lapis::Point(thisTaos.x(k), thisTaos.y(k), thisMhm.crs()))) {
                                    taos.taoVector.addFeature(thisTaos.taoVector.getFeature(k));
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
    return((int)rxs[idx].taos.size());
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
    rxtools::TaoGettersPt g(
        lapis::lico::alwaysAdd<lapis::VectorDataset<lapis::Point>>,
        [](const lapis::ConstFeature<lapis::Point>& ft)->lapis::CoordXY {
            return { ft.getNumericField<lapis::coord_t>("X"), ft.getNumericField<lapis::coord_t>("Y") };
        },
        [](const lapis::ConstFeature<lapis::Point>& ft)->lapis::coord_t {
            return ft.getNumericField<lapis::coord_t>("Height");
        },
        lapis::lico::FixedRadius<lapis::VectorDataset<lapis::Point>>(3),
        [](const lapis::ConstFeature<lapis::Point>& ft)->lapis::coord_t {
            return ft.getNumericField<lapis::coord_t>("Area");
        },
        [dm = dbhModel](const lapis::ConstFeature<lapis::Point>& ft)->double {
            return dm.predict(ft.getNumericField<lapis::coord_t>("Height"),
                lapis::linearUnitPresets::meter, rxtools::linearUnitPresets::centimeter);
        }
    );
    lapis::VectorDataset<lapis::Point> pts;
    pts.addNumericField<lapis::coord_t>("X");
    pts.addNumericField<lapis::coord_t>("Y");
    pts.addNumericField<lapis::coord_t>("Height");
    pts.addNumericField<lapis::coord_t>("Area");

    for (size_t i = 0; i < size / 5; ++i) {
        lapis::Point p = lapis::Point(taoData[5 * i], taoData[5 * i + 1]);
        pts.addGeometry(p);
        pts.back().setNumericField("X", taoData[5 * i]);
        pts.back().setNumericField("Y", taoData[5 * i + 1]);
        pts.back().setNumericField("Height", taoData[5 * i + 3]);
        pts.back().setNumericField("Area", taoData[5 * i + 2]);
    }
    rxs[idx].taos = rxtools::TaoListPt(pts, getters);
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
            outData[i] = (int)rxs[idx].unitMask[i].value();
        }
        else {
            outData[i] = naValue;
        }
    }
};

void setAllometry(double intercept, double slope, int transform) {
    
    dbhModel.parameters.intercept = intercept;
    dbhModel.parameters.slope = slope;
    // 0 = none
    // 1 = square
    // 2 = cube
    // 3 = log-log
    // 4 = suggest
    dbhModel.parameters.transform = (rxtools::allometry::UnivariateLinearModel::Transform)transform;
    dbhModel.inputUnit = lapis::linearUnitPresets::internationalFoot;
    dbhModel.outputUnit = rxtools::linearUnitPresets::inch;

    //the model internally is in  feet and inches, but we want to predict on meters and cm in our getter.
    getters.dbh = [dm = dbhModel, hg = getters.height](const lapis::ConstFeature<lapis::Point>& ft)->double {
        return dm.predict(hg(ft), lapis::linearUnitPresets::meter, rxtools::linearUnitPresets::centimeter);
        };
    allomInit = TRUE;
}

int setAllometryWkt(char* wkt, char* crsWkt, char* fiaPath) {
    auto reader = rxtools::allometry::FIAReader(fiaPath);

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
        std::abort();
    }

    auto dist = lidarDataset->units().value().convertOneToThis(10000, lapis::linearUnitPresets::meter);
    lapis::Extent e(
        poly.boundingBox().xmin() - dist,
        poly.boundingBox().xmax() + dist,
        poly.boundingBox().ymin() - dist,
        poly.boundingBox().ymax() + dist,
        poly.crs());
    std::cout << poly.crs().getPrettyWKT() << "\n";
    auto n = reader.limitByExtent(e);
    if (!n) {
        std::cout << "no fia plots in buffered extent of project area; aborting";
        std::abort();
    }

    reader.makePlotTreeMap(std::vector<std::string>{ "DIA" });
    auto allTrees = reader.collapsePlotTreeMap();
    dbhModel = rxtools::allometry::UnivariateLinearModel(allTrees, "DIA", rxtools::linearUnitPresets::inch);
    getters.dbh = [dm = dbhModel, hg = lidarDataset->heightGetter()](const lapis::ConstFeature<lapis::Point>& ft)->double {
        return dm.predict(hg(ft), lapis::linearUnitPresets::meter, rxtools::linearUnitPresets::centimeter);
        };
    allomInit = TRUE;
    return((int)n);
}

void getAllometry(double* intercept, double* slope, int* transform) {
    *intercept = dbhModel.parameters.intercept;
    *slope = dbhModel.parameters.slope;
    *transform = (int)dbhModel.parameters.transform;
}

//metric in, metric out
void getDbhFromHeight(double* height, double* dbh, int size) {
    for (int i = 0; i < size; ++i) {
        dbh[i] = dbhModel.predict(height[i], lapis::linearUnitPresets::meter, rxtools::linearUnitPresets::centimeter);
    }
}

void getCurrentStructure(int idx, double* ba, double* tph, double* mcs, double* cc) {
    *ba = rxs[idx].currentStructure.ba;
    *tph = rxs[idx].currentStructure.tph;
    *mcs = rxs[idx].currentStructure.mcs;
    *cc = rxs[idx].currentStructure.cc;
}

void calcCurrentStructure(int idx) {
    rxs[idx].currentStructure = rxtools::StructureSummary(
        rxs[idx].taos,
        lapis::Alignment((lapis::Extent)rxs[idx].unitMask, 1, 1),
        rxs[idx].areaHa);
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
    lapis::lico::GraphLico g{ lapis::Alignment((lapis::Extent)rxs[idx].unitMask, 1, 1) };
    g.addDataset(rxs[idx].taos.taoVector, rxs[idx].taos.nodeFactory, lapis::lico::NodeStatus::on);

    for (size_t i = 0; i < g.nodes.size(); ++i) {
        ids[i] = (int)g.nodes[i].findAncestor().index;
    }
}

void makeClumpMap(int idx, int* groupsizes, int* outData, int naValue) {
    try {
        auto b = rxs[idx].basinMap;

        std::unordered_map<int, int> taoIds;

        for (size_t i = 0; i < rxs[idx].taos.size(); ++i) {
            auto e = b.extract(rxs[idx].taos.x(i), rxs[idx].taos.y(i), lapis::ExtractMethod::near);
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
        std::abort();
    }
}

void getSimulatedStructures(int idx, double bbDbh, double* out) {
    try {
        auto rx = rxs[idx];
        auto align = lapis::Alignment((lapis::Extent)rx.unitMask, 1, 1);

        std::vector<rxtools::StructureSummary> structures;

        std::deque<size_t> notBBbase;
        rxtools::TaoListPt base(rx.taos, true);
        for (size_t i = 0; i < rx.taos.size(); i++) {
            if (rx.taos.dbh(i) >= bbDbh) {
                base.taoVector.addFeature(rx.taos.taoVector.getFeature(i));
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
        size_t step = (rx.taos.size() - testTaos.size()) / 10;

        //short to tall.
        for (int i = 1; i < 11; ++i) {
            auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
            char buf[26];
            ctime_s(buf, sizeof(buf), &t);
            std::cout << i << " " << buf << "\n";
            size_t limit = std::min<size_t>(rx.taos.size(), testTaos.size() + step);
            while (testTaos.size() < limit && notBBidx.size()) {
                auto idx = notBBidx.back();
                notBBidx.pop_back();
                testTaos.taoVector.addFeature(rx.taos.taoVector.getFeature(idx));
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
            size_t limit = std::min<size_t>(rx.taos.size(), testTaos.size() + step);
            while (testTaos.size() < limit && notBBidx.size()) {
                auto idx = notBBidx.front();
                notBBidx.pop_front();
                testTaos.taoVector.addFeature(rx.taos.taoVector.getFeature(idx));
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
            size_t limit = std::min<size_t>(rx.taos.size(), testTaos.size() + step);
            while (testTaos.size() < limit && notBBidx.size()) {
                auto idx = notBBidx.back();
                notBBidx.pop_back();
                testTaos.taoVector.addFeature(rx.taos.taoVector.getFeature(idx));
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
        std::abort();
    }
}

void doTreatment(int idx, double dbhMin, double dbhMax) {
    if (!rxInit) {
        std::cerr << "Rx not init\n";
        std::abort();
    }
    try {
        auto trt =  treater.doTreatment(rxs[idx], dbhMin, dbhMax, 10, false, "E:/dropbox/rxgaming paper/treatmentvisual/data/");
        rxs[idx].treatedTaos = std::get<0>(trt);
        rxs[idx].cutTaos = std::get<1>(trt);
        rxs[idx].result = std::get<2>(trt);
        rxs[idx].treatedStructure = rxtools::StructureSummary(
            rxs[idx].treatedTaos,
            lapis::Alignment((lapis::Extent)rxs[idx].unitMask, 1, 1),
            rxs[idx].areaHa);
    }
    catch (std::exception e) {
        std::cout << e.what();
        std::abort();
    }
}

void getTreatedChm(int idx, double* outData, double naValue) {
    auto rx = rxs[idx];
    std::unordered_set<int> basinIds;
    for (size_t i = 0; i < rx.treatedTaos.size(); ++i) {
        auto v = rx.basinMap.extract(rx.treatedTaos.x(i), rx.treatedTaos.y(i), lapis::ExtractMethod::near);
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
    return (int)rxs[idx].treatedTaos.size();
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
    return (int)rxs[idx].cutTaos.size();
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
    lapis::lico::GraphLico g{ lapis::Alignment((lapis::Extent)rxs[idx].unitMask, 1, 1) };
    g.addDataset(rxs[idx].treatedTaos.taoVector, rxs[idx].treatedTaos.nodeFactory, lapis::lico::NodeStatus::on);

    for (size_t i = 0; i < g.nodes.size(); ++i) {
        ids[i] = (int)g.nodes[i].findAncestor().index;
    }
}

void getTreatedClumpMap(int idx, int* inData, int inNaValue, int* groupsizes, int* outData, int naValue) {
    try {
        auto rx = rxs[idx];
        std::unordered_set<int> basinIds;
        for (size_t i = 0; i < rx.treatedTaos.size(); ++i) {
            auto v = rx.basinMap.extract(rx.treatedTaos.x(i), rx.treatedTaos.y(i), lapis::ExtractMethod::near);
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
        std::abort();
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