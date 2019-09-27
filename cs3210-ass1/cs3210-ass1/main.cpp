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
#include "vector2.h"
#include "geometry.h"
#include "Particle.h"

vector2 gStageSize;
double gStepSize;
int gNumSteps;
int gStepNumber = 0;
bool gPrintAll = false;

inline double fRand(double fMin, double fMax)
{
    return fMin + ((double)rand() / RAND_MAX) * (fMax - fMin);
}

inline void PrintParticle(const Particle particle)
{
    std::cout << particle.index << ' ';
    std::cout << gStepNumber << ' ';
    std::cout << particle.position.x << ' ';
    std::cout << particle.position.y << ' ';
    std::cout << particle.velocity.x << ' ';
    std::cout << particle.velocity.y << '\n';
}

int main(int argc, char *argv[])
{
    // Num of particles, Size of square, Radius of particle, and number of steps
    int N, L; double r;

    std::cin >> N >> L >> r >> gNumSteps;
    std::vector<Particle> particles;
    particles.reserve(N);

    gStageSize = vector2(L, L);

    std::string mode;
    std::cin >> mode;
    if (mode == "print") {
        // print for every timestep
        gPrintAll = true;
    }

    // Generate particles based on data
    {
        int particleIndex;
        while (std::cin >> particleIndex) {
            vector2 initialPosition, initialVelocity;
            std::cin >> initialPosition.x >> initialPosition.y >> initialVelocity.x >> initialVelocity.y;
            particles.emplace_back(Particle(initialPosition, initialVelocity, r, particleIndex));
        }
    }

    // Generate random particles
    srand(time(NULL));
    double minVelocity = L / 4;
    double maxVelocity = L / (8 * r);
    while (particles.size() < N) {
        int sign = (rand() % 2) ? 1 : -1;
        vector2 initialPosition(fRand(0.0, L), fRand(0.0, L));
        vector2 initialVelocity(sign * fRand(minVelocity, maxVelocity), sign * fRand(minVelocity, maxVelocity));

        particles.push_back(Particle(initialPosition, initialVelocity, r, particles.size()));
    }

    return 0;
}
