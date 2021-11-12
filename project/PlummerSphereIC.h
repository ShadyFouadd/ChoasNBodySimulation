//
// Created by dia2s on 8/19/2021.

#pragma once

#include "fstream"
#include <iomanip>
#include "M:\Programming Libraries\C++\Boost 66\boost_1_66_0\boost\multi_array.hpp"

//using boost::extents;
typedef boost::multi_array<long double, 2 > Two_D_Array;
void PlumIC(const int seed, const int Num, const long double ZMass, const long double Pmass, long double& TKE,
    Two_D_Array& PositionsL, Two_D_Array& VeL, Two_D_Array& PositionsS, Two_D_Array& VeS);



