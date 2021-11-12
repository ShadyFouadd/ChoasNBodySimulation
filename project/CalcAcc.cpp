//
// Created by dia2s on 8/19/2021.
//

#include "CalcAcc.h"
void CalcAcc(const int Num, const long double pMass, const Two_D_Array& Pos, Two_D_Array& Acc, long double TotPE, One_D_Array PotEn)
{
    long double dS[3]; //dS = {dx, dy, dz}
    long double dS2; //dS^2
    const long double Soft = 0.0005;
    long double Dist; // SQRT(dS^2 + Softening potential)
    long double InvDist; //Dist^(-1)
    long double Coeff; //pMass*Dist^3

    TotPE = 0;
    for (int i = 0; i < Num; ++i)
    {
        for (int j = 0; j < 3; ++j) {
            Acc[i][j] = 0;
        }
        PotEn[i] = 0;
    }

    for (int i = 0; i < Num - 1; ++i)
    {
        for (int j = i + 1; j < Num; ++j)
        {

            dS2 = 0;
            for (int k = 0; k < 3; ++k)
            {
                dS[k] = Pos[i][k] - Pos[j][k];
                dS2 += dS[k] * dS[k];
            }
            Dist = (long double)sqrt(dS2 + Soft);
            InvDist = pow(Dist, -1);
            Coeff = pMass * pow(Dist, -3);

            TotPE -= pow(pMass, 2) * InvDist; //Increments total potential energy by this pair's PE

            for (int k = 0; k < 3; ++k)
            {
                Acc[j][k] += Coeff * dS[k];
                Acc[i][k] -= Coeff * dS[k];
            }
            PotEn[i] -= pow(pMass, 2) * InvDist;
            PotEn[j] -= pow(pMass, 2) * InvDist;
        }
    }
}
