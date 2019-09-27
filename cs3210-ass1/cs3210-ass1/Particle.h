#pragma once

#include "vector2.h"

// constant radius for all particles
double gParticleRadius;

struct Particle
{
    Particle(vector2 pos, vector2 vel, uint32_t index)
        : position(pos), velocity(vel), index(index), numWallCollisions(0), numParticleCollisions(0)
    {
    };

    vector2 position;
    vector2 velocity;

    uint32_t index;

    uint32_t numWallCollisions;
    uint32_t numParticleCollisions;
};
