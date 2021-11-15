//
// Created by dia2s on 8/19/2021.
//

#ifndef MAIN_CPP_CALCACC_H
#define MAIN_CPP_CALCACC_H
#include "boost\multi_array.hpp"

typedef boost::multi_array<long double, 1> One_D_Array;
typedef boost::multi_array<long double, 2> Two_D_Array;

void CalcAcc(const int Num, const long double pMass, const Two_D_Array& Pos, Two_D_Array& Acc, long double TotPE, One_D_Array PotEn);

#endif //MAIN_CPP_CALCACC_H
