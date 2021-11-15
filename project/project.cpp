#include <iostream>
#include "vector"
#include "boost/multi_array.hpp"
#include <time.h>
//#include "fstream"
//#include <iomanip>
#include "PlummerSphereIC.h"
#include "CalcAcc.h"
#include "Analyser.h"
#include "NewStepper.h"
//#include <boost/compute/utility/extents.hpp>

/*
* Why do keep send the same Total Ke and Total PE to stepper? Do we have to make two versions of the total energy? (To Diaa)
* Also, do we need to send the energy by reference?
 *
 * Yes, we need to send them by reference so they get updated. I actually forgot to
 * do that for TotPE and PotEn in CalcAcc too
*/

using namespace std;
using boost::extents;
typedef boost::multi_array<long double, 2 > Two_D_Array;
typedef boost::multi_array<long double, 1 > One_D_Array;
long double pi = (long double)4.0 * atan(1.0);

int main() {
    /////DIAA
    int NLeap = 2;
    int CodeTime = time(0);
    const int Num = 1000;
    int seed = 0;
    long double ZMass = 1.0;
    AnalyserValues ERROR_VALS;
    Two_D_Array PositionL(extents[Num][3]), PositionS(extents[Num][3]);//position vectors of particles
    Two_D_Array VeL(extents[Num][3]), VeS(extents[Num][3]);//velocity vectors of particles (Large ts vs small ts)
    Two_D_Array AccL(extents[Num][3]), AccS(extents[Num][3]);//acceleration vectors of particles
    One_D_Array PEL(extents[Num]), PES(extents[Num]);
    long double TotalPotentialEnergy = 0, TotalKineticEnergy = 0;
    long double TE = TotalPotentialEnergy + TotalKineticEnergy; // TE definition and declaration
    long double pMass = (long double)ZMass / Num;
    long double Tend = 100.0;
    long double dT = 1.0;
    long double DTLARGE = 1.0; //I don't know what is the purpose of DTLARGE in the original code
    //to repeat the simulation for dT and DtLarge (one should be double the other) then compare the error difference

    PlumIC(seed, Num, ZMass, pMass, TotalKineticEnergy, PositionL, VeL, PositionS, VeS); // To Diaa: Where is this in the code?
    //InitialPlummerSphere if I recall correctly. It was the final subroutine in the code

   // CalcAcc(Num, pMass, PositionL, AccL, TotalPotentialEnergy, PEL);
   // ERROR_VALS.Analyser(Num, pMass, PositionL, PositionS, VeL, VeS, PEL, PES);
    cout << time(0) - CodeTime;

    /*
    * Output plummer sphere into a file CHECKED
    * A loop to initialize the system CHECKED
    * Initialise accellerations, energies etc for both small step and big step
    * Initialise differences  (between small step and big step) (Call analyser)
    * Initialise counter (for output) and start looping until we reach the max integration time Tend
    * Do the output
    */
    double T = 0;
//    ofstream myFile;
//    myFile.open("Plummer Output.txt");
//
//    for (int i = 0; i < Num; i++)
//    {
//        cout << i + 1 << "." << PositionS[i][1] << "," << PositionS[i][2] << "," << PositionS[i][3] << endl;
//    }
//    myFile.close();
// this process is already executed inside ICPlum

    //A loop to initialize the system 

    for (int i = 0; i < Num; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            PositionL[i][j] = PositionS[i][j];
            VeL[i][j] = VeS[i][j];
            TE += 0.5 * pMass * VeS[i][j] * VeS[i][j];
        }
    }

    //Initialise accellerations, energies etc for both small step and big step

    CalcAcc(Num, pMass, PositionS, AccS, TotalPotentialEnergy, PES);
    CalcAcc(Num, pMass, PositionL, AccL, TotalPotentialEnergy, PEL);

    //Initialise differences  (between small step and big step) (Call analyser)

    ERROR_VALS.Analyser(Num, pMass, PositionL, PositionS, VeL, VeS, PEL, PES);

    //Here the original code defines ENL as TE plus PE but I don't get that. It also defines the virial ration, what the fuck is that?

    //Initialise counter(for output) and start looping until we reach the max integration time Tend
    int Counter = 0;
    while (T < Tend)
    {
        //For the larger step
        dT = DTLARGE;
        double dT2 = dT * dT;
        PlumStepper(dT, dT2, Num, pMass, PositionL, VeL, TotalPotentialEnergy, TotalKineticEnergy, PEL, AccL);


        // And now for smaller time step
        dT = DTLARGE / NLeap; // I don't know why we divide by NLeap here
        dT2 = dT * dT;

        for (int i = 0; i < NLeap; i++)
            PlumStepper(dT, dT2, Num, pMass, PositionS, VeS, TotalPotentialEnergy, TotalKineticEnergy, PES, AccS);
        T += DTLARGE;
        Counter += 1;

        if (Counter % 10 == 0) //Why?
            ERROR_VALS.Analyser(Num, pMass, PositionL, PositionS, VeL, VeS, PEL, PES);

    }
}
