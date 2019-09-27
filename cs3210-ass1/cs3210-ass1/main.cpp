/**
    CS3210 Assignment 1 Part 1 (OpenMP)
    Authors: Kuan Wei Heng (A0121712X)
             Chong Jun Hong, Dominic (A0121237U)
    std=C++11
*/
#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "point.h"
#include "geometry.h"
#include "Particle.h"

point gStageSize;
double gStepSize;
int gNumSteps;
int gStepNumber = 0;
bool gPrintAll = false;

std::vector<Particle> particles;

inline double fRand(double fMin, double fMax)
{
    return fMin + ((double)rand() / RAND_MAX) * (fMax - fMin);
}

inline void PrintParticle(const Particle particle);

inline void PrintParticle(const Particle particle)
{
    printf("%i %i %10.8f %10.8f %10.8f\n %10.8f\n\n",
           particle.index,
           gStepNumber,
           particle.position.x,
           particle.position.y,
           particle.velocity.x,
           particle.velocity.y
    );
}

int main(int argc, char *argv[])
{
    // Num of particles, Size of square, Radius of particle, and number of steps
    int N, L;
    scanf_s("%i", &N);
    scanf_s("%i", &L);
    gStageSize = point(L, L);
    scanf_s("%lf", &gParticleRadius);
    scanf_s("%i", &gNumSteps);

    std::string inputBuffer;
    scanf_s("%s", inputBuffer);
    if (inputBuffer.compare("print")) {
        // print for every timestep
        gPrintAll = true;
    }

    // Generate particles based on data
    int particleCount = 0;
    int particleIndex = 0;
    particles.reserve(N);

    point initialPosition, initialVelocity;
    while (scanf_s("%i %lf %lf %lf %lf", &particleIndex,
           &initialPosition.x,
           &initialPosition.y,
           &initialVelocity.x,
           &initialVelocity.y) == 5) {
        // TODO: Create particle
        particles.push_back(Particle(initialPosition, initialVelocity, particleIndex));
        particleCount++;
    }

    // Generate random particles
    srand(time(NULL));
    double minVelocity = L / 4;
    double maxVelocity = L / (8 * gParticleRadius);
    for (; particleCount < N; particleCount++) {
        // TODO: Create random particles
        int sign = (rand() % 2) ? 1 : -1;
        initialPosition = point(fRand(0.0, L), fRand(0.0, L));
        initialVelocity = point(sign * fRand(minVelocity, maxVelocity), sign * fRand(minVelocity, maxVelocity));

        particles.push_back(Particle(initialPosition, initialVelocity, particleCount));
    }

    return 0;
}
