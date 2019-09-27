#pragma once

#include "point.h"

// constant radius for all particles
double gParticleRadius;

struct Particle
{
    Particle(point pos, point vel, uint32_t index)
        : position(pos), velocity(vel), index(index), numWallCollisions(0), numParticleCollisions(0)
    {
    };

    point position;
    point velocity;

    uint32_t index;

    uint32_t numWallCollisions;
    uint32_t numParticleCollisions;
};
