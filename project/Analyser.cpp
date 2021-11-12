//
// Created by dia2s on 8/19/2021.
//
#include <math.h>
#include "Analyser.h"
void ErrorValues::Regression()
{
    rErr = sqrt(error / eNorm);
}
void ErrorValues::AnalyseVector(const int Num, const Two_D_Array& Small, const Two_D_Array& Large)
{
    for (int i = 0; i < Num; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            this->error += pow((Small[i][j] - Large[i][j]), 2);
            this->eNorm += pow(Small[i][j], 2);
        }
    }
    this->Regression();
}
void ErrorValues::AnalysePotential(const int Num, long double pMass, const Two_D_Array& PosS,
    const Two_D_Array& PosL, const Two_D_Array& VeS, const Two_D_Array& VeL)
{
    long double PXS, PYS, PZS, PXL, PYL, PZL;
    for (int i = 0; i < Num; ++i)
    {
        PXS = pMass * (PosS[i][1] * VeS[i][2] - PosS[i][2] * VeS[i][1]);
        PYS = pMass * (PosS[i][2] * VeS[i][0] - PosS[i][0] * VeS[i][2]);
        PZS = pMass * (PosS[i][0] * VeS[i][1] - PosS[i][1] * VeS[i][0]);
        PXL = pMass * (PosL[i][1] * VeL[i][2] - PosL[i][2] * VeL[i][1]);
        PYL = pMass * (PosL[i][2] * VeL[i][0] - PosL[i][0] * VeL[i][2]);
        PZL = pMass * (PosL[i][0] * VeL[i][1] - PosL[i][1] * VeL[i][0]);

        this->error += pow((PXS - PXL), 2) +
            pow((PYS - PYL), 2) + pow((PZS - PZL), 2);
        this->eNorm += pow(PXS, 2) + pow(PYS, 2) + pow(PZS, 2);
    }
    this->Regression();
}
void ErrorValues::AnalyseTotalEnergy(const int Num, long double pMass, const One_D_Array& PES,
    const One_D_Array& PEL, const Two_D_Array& VeS, const Two_D_Array& VeL)
{
    long double EnergyS, EnergyL;
    for (int i = 0; i < Num; ++i)
    {
        EnergyS = PES[i];
        EnergyL = PEL[i];
        for (int j = 0; j < 3; ++j)
        {
            EnergyS += 0.5 * pMass * pow(VeS[i][1], 2);
            EnergyL += 0.5 * pMass * pow(VeL[i][1], 2);
        }
        this->error += pow((EnergyL - EnergyS), 2);
        this->eNorm += pow(EnergyS, 2);
    }
    this->Regression();
}

void AnalyserValues::Analyser(const int Num, long double pMass, const Two_D_Array& PosL, const Two_D_Array& PosS,
    const Two_D_Array& VeL, const Two_D_Array& VeS, const One_D_Array& PEL, const One_D_Array& PES)
{
    C.AnalyseVector(Num, PosS, PosL);
    V.AnalyseVector(Num, VeS, VeL);
    P.AnalysePotential(Num, pMass, PosL, PosS, VeL, VeS);
    E.AnalyseTotalEnergy(Num, pMass, PES, PEL, VeL, VeS);
}