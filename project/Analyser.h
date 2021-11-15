#ifndef MAIN_CPP_ANALYSER_H
#define MAIN_CPP_ANALYSER_H
#include "boost\multi_array.hpp"


typedef boost::multi_array<long double, 1> One_D_Array;
typedef boost::multi_array<long double, 2> Two_D_Array;


struct ErrorValues
{
    long double error = 0, eNorm = 0, rErr = 0;
    void Regression();
    void AnalyseVector(int Num, const Two_D_Array& Small, const Two_D_Array& Large);
    void AnalysePotential(int Num, long double pMass, const Two_D_Array& PosS,
        const Two_D_Array& PosL, const Two_D_Array& VeS, const Two_D_Array& VeL);
    void AnalyseTotalEnergy(int Num, long double pMass, const One_D_Array& PES,
        const One_D_Array& PEL, const Two_D_Array& VeS, const Two_D_Array& VeL);
};

struct AnalyserValues
{
    ErrorValues C, V, P, E; //naming should be adjusted
    void Analyser(int Num, long double pMass, const Two_D_Array& PosL, const Two_D_Array& PosS,
        const Two_D_Array& VeL, const Two_D_Array& VeS, const One_D_Array& PEL, const One_D_Array& PES);
};

#endif //MAIN_CPP_ANALYSER_H
#pragma once
