/*
Copyright(C) 2024  University of Washington
This program is free software : you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.If not, see https ://www.gnu.org/licenses/.

Bryce Bartl - Geller
University of Washington Forest Resilience Lab
12 / 6 / 2024

cpprxgaming.h
This is a .dll that does some of the slower things from the rxgaming tool faster here.
This is mainly running the treatment simulation and loading and manipulating data from the larger
    lidar datsets.
*/

#pragma once

#ifdef CPPRXGAMING_EXPORTS
#define CPPRXGAMING_API __declspec(dllexport)
#else
#define CPPRXGAMING_API __declspec(dllimport)
#endif

#include<string>
#include<mutex>
#include <unordered_set>
#include "raster.hpp"
#include "graphlico.hpp"
#include "treatment.hpp"
#include "rxunit.hpp"


class RxGamingRxUnit : public rxtools::RxUnit {
public:
    lapis::Raster<int> mhm;
    lapis::Raster<double> chm;
    lapis::Raster<int> basinMap;
    rxtools::TaoListPt cutTaos;
    rxtools::treatmentResult result;

    RxGamingRxUnit(lapis::Raster<lapis::cell_t> mask, lapis::Raster<int> mhm, lapis::Raster<double> chm, lapis::Raster<int> basinMap,
        rxtools::TaoListPt taos, double convFactor) :  RxUnit(mask, taos) {
        this->mhm = mhm;
        this->chm = chm;
        this->basinMap = basinMap;
    };

    RxGamingRxUnit() = default;
};

extern "C" CPPRXGAMING_API void setProjDataDirectory(const char* searchpath);

extern "C" CPPRXGAMING_API void setSeed(int n);

extern "C" CPPRXGAMING_API bool initLidarDataset(char* path);

extern "C" CPPRXGAMING_API void reprojectPolygon(char* wkt, char* crsWkt, char* outWkt);

extern "C" CPPRXGAMING_API bool queueRx(char* wkt, char* crsWkt);

extern "C" CPPRXGAMING_API void doPreProcessing(int nThread);
void rastersAndTaosThread(const int nThread, const int thisThread, std::mutex& mut, int& sofar, const double canopycutoff,
    const double coregapdist, std::pair<lapis::coord_t, lapis::coord_t> expectedRes);

extern "C" CPPRXGAMING_API int nTaos(int idx);
extern "C" CPPRXGAMING_API void getTaos(int idx, double* outData);
extern "C" CPPRXGAMING_API void setTaos(int idx, double* taoData, int size);

extern "C" CPPRXGAMING_API void setMhm(int idx, int* data, int nrow, int ncol, double xres, double yres, double xmin, double ymin, int naValue, char* projWKT);
extern "C" CPPRXGAMING_API void getMhmMeta(int idx, int* nrow, int* ncol, double* xres, double* yres, double* xmin, double* ymin, char* projWKT);
extern "C" CPPRXGAMING_API void getMhm(int idx, int* outData, int naValue);

extern "C" CPPRXGAMING_API void setChm(int idx, double* data, int nrow, int ncol, double xres, double yres, double xmin, double ymin, double naValue, char* projWKT);
extern "C" CPPRXGAMING_API void getChmMeta(int idx, int* nrow, int* ncol, double* xres, double* yres, double* xmin, double* ymin, char* projWKT);
extern "C" CPPRXGAMING_API void getChm(int idx, double* outData, double naValue);

extern "C" CPPRXGAMING_API void setBasin(int idx, int* data, int nrow, int ncol, double xres, double yres, double xmin, double ymin, int naValue, char* projWKT);
extern "C" CPPRXGAMING_API void getBasinMeta(int idx, int* nrow, int* ncol, double* xres, double* yres, double* xmin, double* ymin, char* projWKT);
extern "C" CPPRXGAMING_API void getBasin(int idx, int* outData, int naValue);

extern "C" CPPRXGAMING_API void setMask(int idx, int* data, int nrow, int ncol, double xres, double yres, double xmin, double ymin, int naValue, char* projWKT);
extern "C" CPPRXGAMING_API void getMaskMeta(int idx, int* nrow, int* ncol, double* xres, double* yres, double* xmin, double* ymin, char* projWKT);
extern "C" CPPRXGAMING_API void getMask(int idx, int* outData, int naValue);

extern "C" CPPRXGAMING_API void setAllometryFiaPath(char* path);
extern "C" CPPRXGAMING_API void setAllometry(double intercept, double slope, int transform);
extern "C" CPPRXGAMING_API int setAllometryWkt(char* wkt, char* crsWkt);
extern "C" CPPRXGAMING_API void getAllometry(double* intercept, double* slope, int* transform);
extern "C" CPPRXGAMING_API void getDbhFromHeight(double* height, double* dbh, int size);

extern "C" CPPRXGAMING_API void getCurrentStructure(int idx, double* ba, double* tph, double* mcs, double* osi, double* cc);
extern "C" CPPRXGAMING_API void calcCurrentStructure(int idx);
extern "C" CPPRXGAMING_API void setTargetStructure(int idx, double ba, double tph, double mcs, double osi, double cc);
extern "C" CPPRXGAMING_API void getTargetStructure(int idx, double* ba, double* tph, double* mcs, double* osi, double* cc);

extern "C" CPPRXGAMING_API void getRawClumps(int idx, int* ids);
extern "C" CPPRXGAMING_API void makeClumpMap(int idx, int* groupsizes, int* outData, int naValue);

extern "C" CPPRXGAMING_API void getSimulatedStructures(int idx, double bbDbh, double* out);

//Do treatment on currently loaded Taos
extern "C" CPPRXGAMING_API void doTreatment(int idx, double dbhMin, double dbhMax);
extern "C" CPPRXGAMING_API void getTreatedChm(int idx, double* outData, double naValue);
extern "C" CPPRXGAMING_API void getTreatedBasin(int idx, int* outData, int naValue);
extern "C" CPPRXGAMING_API int getNTreatedTaos(int idx);
extern "C" CPPRXGAMING_API void getTreatedTaos(int idx, double* outData);
extern "C" CPPRXGAMING_API int getNCutTaos(int idx);
extern "C" CPPRXGAMING_API void getCutTaos(int idx, double* outData);
extern "C" CPPRXGAMING_API int getTreatmentResult(int idx);
extern "C" CPPRXGAMING_API void getTreatedStructure(int idx, double* ba, double* tph, double* mcs, double* osi, double* cc);
extern "C" CPPRXGAMING_API void getTreatedRawClumps(int idx, int* ids);
extern "C" CPPRXGAMING_API void getTreatedClumpMap(int idx, int* inData, int inNaValue, int* groupsizes, int* out, int naValue);

//export info needed to reconstruct without preprocessing.
extern "C" CPPRXGAMING_API void setRxsSize(int i);

extern "C" CPPRXGAMING_API void setConvFactor(double cf);
extern "C" CPPRXGAMING_API void getConvFactor(double* cf);