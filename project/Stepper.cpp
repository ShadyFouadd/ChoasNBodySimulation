#include "Stepper.h"

void STEPPER(double DT, double DT2, const int Num, const long double pMass, Two_D_Array& X, Two_D_Array& V,
	long double TotPE, long double TotKE, One_D_Array PotEn, const Two_D_Array& A)
{
	//double TotKE = 0;
	Two_D_Array AOLD(extents[Num][3]);

	for (int i = 0; i < Num; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			X[i][j] += DT * V[i][j] + 0.5 * DT2 * A[i][j];
			// TotKE += 0.5 * pMass * V[i][j] * V[i][j]; //This was commented in the original code.
			AOLD[i][j] = A[i][j];
			// CalcAcc(Num, pMass, X, A, TotPE, PotEn);
		}
	}
//	CalcAcc(Num, pMass, X, A, TotPE, PotEn);

	for (int i = 0; i < Num; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			V[i][j] += 0.5 * DT * (A[i][j] + AOLD[i][j]); //I should add AOld here 
			TotKE = 0.5 * pMass * V[i][j] * V[i][j];
		}
	}
}