//
// Created by dia2s on 8/19/2021.
//

#include "PlummerSphereIC.h"
static long double pi = (long double)4.0 * atan(1.0);

void PlumIC(const int seed, const int Num, const long double ZMass, const long double Pmass, long double& TKE, Two_D_Array& PositionsL, Two_D_Array& VeL, Two_D_Array& PositionsS, Two_D_Array& VeS)
{
    long double Rmass;
    long double radius = 1000;
    long double velocity;
    long double position[3]; //resembles x,y,z
    long double velocityVec[3];
    long double CM[3] = { 0,0,0 };
    long double CMTrajectory[3] = { 0,0,0 };
    long double phi;
    long double theta;
    long double R1;
    long double R2;
    srand(seed);
    std::ofstream myfile;
    myfile.open("PlummerSphereInitial.txt");
    srand(seed);
    int linew1 = 10, linew2 = 10;
    for (int i = 0; i < Num; ++i)
    {
        radius = 1000;
        while (radius > 10.0)
        {
            Rmass = (long double)rand() / RAND_MAX;
            while (Rmass == 1)
                Rmass = (long double)rand() / RAND_MAX;
            radius = pow((pow(Rmass, -0.6666667) - 1.0), -0.5);
        }
        theta = (long double)pi * rand() / RAND_MAX;
        phi = (long double)2.0 * pi * rand() / RAND_MAX;
        position[0] = radius * cos(phi) * sin(theta);
        position[1] = radius * sin(phi) * sin(theta);
        position[2] = radius * cos(theta);

        //        now we accquire the velocity component
        R1 = (long double)rand() / RAND_MAX;
        R2 = (long double)rand() / RAND_MAX;
        while (pow(R1, 2.0) * pow(1.0 - pow(R1, 2.0), 3.5) < 0.1 * R2)
        {
            R1 = (long double)rand() / RAND_MAX;
            R2 = (long double)rand() / RAND_MAX;
        }
        velocity = R1 * (pow(2.0, 0.5) * pow(1 + pow(radius, 2.0), -0.25));
        theta = (long double)pi * rand() / RAND_MAX;
        phi = (long double)2.0 * pi * rand() / RAND_MAX;
        velocityVec[0] = velocity * cos(phi) * sin(theta);
        velocityVec[1] = velocity * sin(phi) * sin(theta);
        velocityVec[2] = velocity * cos(theta);

        for (int j = 0; j < 3; ++j)
        {
            PositionsL[i][j] = position[j];
            PositionsS[i][j] = position[j];
            VeL[i][j] = velocityVec[j];
            VeS[i][j] = velocityVec[j];

            CM[j] += Pmass * position[j];
            CMTrajectory[j] += Pmass * velocityVec[j];
        }
    }
    double SPos = (3.0 * pi) / 16.0;
    double SVel = sqrt(ZMass / SPos); //ZMass= 1
    for (int i = 0; i < Num; ++i) {
        for (int j = 0; j < 3; ++j) {
            PositionsL[i][j] = SPos * (PositionsL[i][j] - CM[j] / ZMass); //ZMass= 1
            PositionsS[i][j] = SPos * (PositionsS[i][j] - CM[j] / ZMass); //ZMass= 1
            VeL[i][j] = SVel * (VeL[i][j] - CMTrajectory[j] / ZMass);
            VeS[i][j] = SVel * (VeS[i][j] - CMTrajectory[j] / ZMass);
            TKE += 0.5 * Pmass * VeS[i][j] * VeS[i][j];
        }
        myfile << PositionsL[i][0] << std::setw(linew1 - 4) << "" << std::setw(linew2)
            << PositionsL[i][1] << std::setw(linew1 - 4) << "" << std::setw(linew2)
            << PositionsL[i][2] << std::setw(linew1 - 4) << "" << std::setw(linew2)
            << std::endl;
    }
    myfile.close();
}