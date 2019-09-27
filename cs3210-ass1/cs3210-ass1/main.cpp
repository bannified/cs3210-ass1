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
#include <map>
#include <algorithm>
#include "vector2.h"
#include "geometry.h"
#include "Particle.h"

vector2 gStageSize;
double gStepSize;
int gNumSteps;
int gStepNumber = 0;
bool gPrintAll = false;

struct Collision
{
    Collision(uint32_t index1, uint32_t index2, double stepValue)
        : index1(index1), index2(index2), stepValue(stepValue) {}

    uint32_t index1;
    uint32_t index2;
    double stepValue;
};

inline void PrintParticle(const Particle particle);

inline double fRand(double fMin, double fMax)
{
    return fMin + ((double)rand() / RAND_MAX) * (fMax - fMin);
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

    std::vector< Collision > collisionResults;
    collisionResults.reserve(N);

    // Start simulation for gNumSteps
    for (; gStepNumber < gNumSteps; gStepNumber++) {
        // Print all particles
        for (const Particle particle : particles) {
            PrintParticle(particle);
        }

        // Checking for particle-to-particle collision
        for (const Particle& particle : particles) {
            for (const Particle& target : particles) {
                double step = canParticlesCollide(particle, target);
                if (step >= 0) {
                    collisionResults.push_back({particle.index, target.index, step});
                }
            }
        }

        // Checking for particle-to-wall collision
        for (const Particle& particle : particles) {
            // TODO: check wall collision for every particle
            /*double collisionCheckResult = canParticlesCollide(particle, target);
            if (collisionCheckResult >= 0) {
                collisionCheckResults.push_back({ particle.index, target.index, collisionCheckResult });
            }*/
        }

        // Sort collision check results
        // TODO: OMP parallelization
        std::sort(collisionResults.begin(), collisionResults.end(), 
        [](const Collision c1, const Collision c2) {
            return c1.stepValue > c2.stepValue;
        });

        // Resolution
        std::vector<bool> resolved;
        resolved.resize(N, false);
        for (const Collision res : collisionResults) {
            if (resolved[res.index1]) {
                continue;
            }

            if (res.index2 < 0) {
                // TODO: resolve particle-to-wall collision
                resolved[res.index1] = true;
            }
            else {
                if (resolved[res.index2]) {
                    continue;
                }

                resolveP2PCollision(particles[res.index1], particles[res.index2], res.stepValue, gStageSize);
                resolved[res.index1] = true;
                resolved[res.index2] = true;
            }
        }
    }

    for (const Particle particle : particles)
    {
        PrintParticle(particle);
    }

    return 0;
}

inline void PrintParticle(const Particle particle)
{
    std::cout << gStepNumber << ' ';
    std::cout << particle.index << ' ';
    std::cout << particle.position.x << ' ';
    std::cout << particle.position.y << ' ';
    std::cout << particle.velocity.x << ' ';
    std::cout << particle.velocity.y << '\n';
}
