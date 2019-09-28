/**
    CS3210 Assignment 1 Part 1 (OpenMP)
    Authors: Kuan Wei Heng (A0121712X)
             Chong Jun Hong, Dominic (A0121237U)
    std=C++11
*/
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <map>
#include <algorithm>
#include "vector2.h"
#include "collision.h"
#include "Particle.h"
#include <omp.h>

int threads;
vector2 gStageSize;
double gStepSize;
int gNumSteps;
int gStepNumber = 0;
bool gPrintAll = false;

inline void PrintParticle(const Particle particle);

inline double fRand(double fMin, double fMax)
{
    return fMin + ((double)rand() / RAND_MAX) * (fMax - fMin);
}

int main(int argc, char *argv[])
{
    std::cout << "Usage: " << argv[0] << " <threads>\n";

    if (argc >= 2)
        threads = atoi(argv[2]);
    else
        threads = -1;

    // Multiply the matrices
    if (threads != -1) {
        omp_set_num_threads(threads);
    }

#pragma omp parallel
    {
        threads = omp_get_num_threads();
    }

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
        vector2 initialPosition(fRand(r, L-r), fRand(r, L-r));
        vector2 initialVelocity(sign * fRand(minVelocity, maxVelocity), sign * fRand(minVelocity, maxVelocity));

        particles.push_back(Particle(initialPosition, initialVelocity, r, particles.size()));
    }

    // Start simulation for gNumSteps
    for (; gStepNumber < gNumSteps; gStepNumber++) {
        std::vector<Collision> collisionResults;
        collisionResults.reserve(N);

        // Print all particles
        for (const Particle particle : particles) {
            PrintParticle(particle);
        }

        // Checking for particle-to-particle collision
        for (int i = 0; i < particles.size(); i++) {
            const Particle& particle = particles[i];
            for (int j = i + 1; j < particles.size(); j++) {
                const Particle& target = particles[i + 1];
                if (&particle == &target) {
                    continue;
                }
                double step = canParticlesCollide(particle, target);
                if (isStepValid(step)) {
                    collisionResults.push_back({particle.index, target.index, step});
                }
            }
        }

        // Checking for particle-to-wall collision
        for (const Particle& particle : particles) {
            Collision result = detectWallCollision(particle, gStageSize);
            if (isStepValid(result.stepValue)) {
                collisionResults.push_back(result);
            }
        }

        // Sort collision check results
        // TODO: OMP parallelization
        std::sort(collisionResults.begin(), collisionResults.end());

        // Resolve collisions starting from earliest
        std::vector<bool> resolved(N);
        for (const Collision res : collisionResults) {
            if (resolved[res.index1]) continue;
            if (res.index2 < 0) {
                resolveWallCollision(particles[res.index1], res.index2, res.stepValue, gStageSize);
                clamp(particles[res.index1], gStageSize);
                resolved[res.index1] = true;
            } else {
                if (resolved[res.index2]) continue;
                resolveP2PCollision(particles[res.index1], particles[res.index2], res.stepValue);
                clamp(particles[res.index1], gStageSize);
                clamp(particles[res.index2], gStageSize);
                resolved[res.index1] = true;
                resolved[res.index2] = true;
            }
        }

        // move remaining particles
        for (int i=0; i < N; i++) {
            if (!resolved[i]) {
                particles[i].position += particles[i].velocity;
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
    std::cout << std::fixed << std::setprecision(8);
    std::cout << gStepNumber << ' ';
    std::cout << particle.index << ' ';
    std::cout << particle.position.x << ' ';
    std::cout << particle.position.y << ' ';
    std::cout << particle.velocity.x << ' ';
    std::cout << particle.velocity.y << '\n';
}
