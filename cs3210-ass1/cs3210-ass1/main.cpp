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

inline void PrintParticle(const Particle p);

inline double fRand(double fMin, double fMax)
{
    return fMin + ((double)rand() / RAND_MAX) * (fMax - fMin);
}

int main(int argc, char *argv[])
{
    std::cout << "Usage: " << argv[0] << " <threads>\n";

    if (argc >= 2)
        threads = atoi(argv[1]);
    else
        threads = -1;

    if (threads != -1) {
        omp_set_num_threads(threads);
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
    while ((int)particles.size() < N) {
        int sign = (rand() % 2) ? 1 : -1;
        vector2 initialPosition(fRand(r, L-r), fRand(r, L-r));
        vector2 initialVelocity(sign * fRand(minVelocity, maxVelocity), sign * fRand(minVelocity, maxVelocity));

        particles.push_back(Particle(initialPosition, initialVelocity, r, particles.size()));
    }

    // Start simulation for gNumSteps
    for (; gStepNumber < gNumSteps; gStepNumber++) {
        std::vector<Collision> collisions;
        collisions.reserve(N);

        // Print all particles
        for (const Particle p : particles) PrintParticle(p);

        // Checking for particle-to-particle collision
        // Each thread works on one particle's checking
#pragma omp parallel for shared(particles, collisions)
        for (size_t i = 0; i < particles.size(); i++) {
            const Particle& particle = particles[i];
            for (size_t j = i + 1; j < particles.size(); j++) {
                const Particle& target = particles[j];
                if (&particle == &target) {
                    continue;
                }
                double step = detectParticleCollision(particle, target);
                if (isStepValid(step)) {
                    collisions.push_back({particle.index, target.index, step});
                }
            }
        }

        // Checking for particle-to-wall collision
        for (const Particle& particle : particles) {
            Collision result = detectWallCollision(particle, gStageSize);
            if (isStepValid(result.stepValue)) {
                collisions.push_back(result);
            }
        }

        // Sort collision check results
        std::sort(collisions.begin(), collisions.end());

        std::vector<bool> resolved(N);

        std::vector<Collision> validatedResults; // Stores all of the collision results to be resolved
        validatedResults.reserve(N / 2);

        // Pick out the collisions that are valid starting from the smallest step value
        // and not allowing for repeated collisions for any one particle.
        for (size_t i = 0; i < collisions.size(); i++) {
            Collision res = collisions[i];
            if (resolved[res.index1]) continue;
            if (res.index2 < 0) {
                validatedResults.push_back(res);
                resolved[res.index1] = true;
            } else {
                if (resolved[res.index2]) continue;
                validatedResults.push_back(res);
                resolved[res.index1] = true;
                resolved[res.index2] = true;
            }
        }

        // Resolve the validated collisions in parallel
#pragma omp parallel for shared(particles, validatedResults)
        for (size_t i = 0; i < validatedResults.size(); i++) {
            Collision res = validatedResults[i];
            if (res.index2 < 0) {
                resolveWallCollision(particles[res.index1], res.index2, res.stepValue, gStageSize);
                clamp(particles[res.index1], gStageSize);
            } else {
                resolveParticleCollision(particles[res.index1], particles[res.index2], res.stepValue);
                clamp(particles[res.index1], gStageSize);
                clamp(particles[res.index2], gStageSize);
            }
        }

        // move remaining particles
#pragma omp parallel for shared(particles, resolved)
        for (int i=0; i < N; i++) {
            if (!resolved[i]) {
                particles[i].position += particles[i].velocity;
            }
        }
    }

    for (const Particle p : particles) PrintParticle(p);

    return 0;
}

inline void PrintParticle(const Particle p)
{
    std::cout << std::fixed << std::setprecision(8);
    std::cout << gStepNumber << ' ';
    std::cout << p.index << ' ';
    std::cout << p.position.x << ' ';
    std::cout << p.position.y << ' ';
    std::cout << p.velocity.x << ' ';
    std::cout << p.velocity.y << '\n';
}
