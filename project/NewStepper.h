//
// Created by dia2s on 11/15/2021.
//
#ifndef CHOASNBODYSIMULATION_NEWSTEPPER_H
#define CHOASNBODYSIMULATION_NEWSTEPPER_H

#include "boost/multi_array.hpp"
#include "CalcAcc.h"

using boost::extents;
typedef boost::multi_array<long double, 1 > One_D_Array;
typedef boost::multi_array<long double, 2 > Two_D_Array;

void PlumStepper(double DT, double DT2, const int Num, const long double pMass, Two_D_Array& X, Two_D_Array& V,
             long double TotPE, long double TotKE, One_D_Array PotEn, Two_D_Array& A);


#endif //CHOASNBODYSIMULATION_NEWSTEPPER_H
