#pragma once

#include "vector2.h"

struct Particle
{
    Particle()
    {
        position = vector2(0, 0);
        velocity = vector2(0, 0);
        radius = 0;
        index = -5;
        numWallCollisions = 0;
        numParticleCollisions = 0;
    }

    Particle(vector2 pos, vector2 vel, double r, uint32_t index)
        : position(pos), velocity(vel), radius(r), index(index), numWallCollisions(0), numParticleCollisions(0)
    {
    };

    vector2 position;
    vector2 velocity;
    double radius;

    uint32_t index;

    uint32_t numWallCollisions;
    uint32_t numParticleCollisions;
};
