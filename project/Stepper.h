#pragma once
//#include "boost/multi_array.hpp"
#include "CalcAcc.h"

typedef boost::multi_array<long double, 1 > One_D_Array;
typedef boost::multi_array<long double, 2 > Two_D_Array;

void STEPPER(double DT, double DT2, const int Num, const long double pMass, Two_D_Array& X, Two_D_Array& V,
	long double TotPE, long double TotKE, One_D_Array PotEn, const Two_D_Array& A);