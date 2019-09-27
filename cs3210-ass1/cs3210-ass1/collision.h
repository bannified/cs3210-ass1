#pragma once

#include "vector2.h"
#include "Particle.h"

struct Collision
{
    Collision(uint32_t index1, uint32_t index2, double stepValue)
        : index1(index1), index2(index2), stepValue(stepValue) {}

    uint32_t index1;
    uint32_t index2;
    double stepValue;

    bool operator<(const Collision& rhs) const
    {
        return stepValue < rhs.stepValue;
    }
};


// find closest point on line segment pq to point r
vector2 closestPointOnLine(vector2 p, vector2 q, vector2 r){
    double x1 = p.x;
    double y1 = p.y;

    double x2 = q.x;
    double y2 = q.y;

    double x0 = r.x;
    double y0 = r.y;

    double a1 = y2 - y1;
    double b1 = x1 - x2;
    double c1 = (y2 - y1) * x1 + (x1 - x2) * y1;
    double c2 = -b1 * x0 + a1 * y0;
    double det = a1 * a1 + b1 * b1;
    if (det != 0) {
        return { (a1 * c1 - b1 * c2) / det, (a1 * c2 + b1 * c1) / det };
    } else {
        return { x0, y0 };
    }
}

/**
 * Checks if two moving particles can collide.
 * If they do, the factor/step (between 0.0 and 1.0) of the timestep is returned. 
 * Otherwise, a negative value is returned.
 */
double canParticlesCollide(const Particle& a, const Particle& b) {
    double distance = dist(b.position, a.position);
    double sumRadii = a.radius + b.radius;
    distance -= sumRadii;

    vector2 resultVector = a.velocity - b.velocity;
    double resultMag = magnitude(resultVector);

    // Early escape: Can't reach just based on maximum travel possible
    if (resultMag < distance) {
        return -1;
    }

    vector2 unitResultVector = resultVector;
    unitResultVector.normalize();

    vector2 c = b.position - a.position;
    double d = unitResultVector * c;

    // Early escape: result vector does not cause the particles to move closer together.
    // since one is not moving in the direction of the other.
    if (d <= 0) {
        return -1;
    }

    double lengthC = magnitude(c);
    double fSquared = (lengthC * lengthC) - d;

    double sumRadiiSquared = sumRadii * sumRadii;

    // Escape: closest that a will get to b.
    if (fSquared >= sumRadiiSquared) {
        return -1;
    }

    double tSquared = sumRadiiSquared - fSquared;

    // negative tSquared. Probably don't have to do this check because the one preceding 
    // this one already ensures that tSquared isn't negative.
    if (tSquared < 0) {
        return -1;
    }

    double distanceToCollide = d - std::sqrt(tSquared);

    // Ensure that distance A has to move to touch B
    // is not greater than the magnitude of the movement vector
    if (resultMag < distanceToCollide) {
        return -1;
    }

    // the final displacement that the particle would have just before the collision.
    // can also choose to return this in a result vector;
    vector2 finalVector = unitResultVector * distanceToCollide;

    return magnitude(finalVector) / resultMag;
}

// keep a particle within bounds
void clamp(Particle& p, vector2 stageSize) {
    p.position.x = std::min(stageSize.x - p.radius, p.position.x);
    p.position.x = std::max(p.radius, p.position.x);
    p.position.y = std::min(stageSize.y - p.radius, p.position.y);
    p.position.y = std::max(p.radius, p.position.y);
}

void resolveP2PCollision(Particle& a, Particle& b, double stepProportion)
{
    vector2 aImpact = a.position + a.velocity * stepProportion;
    vector2 bImpact = b.position + b.velocity * stepProportion;

    double d = std::sqrt(std::pow(aImpact.x + bImpact.x, 2) + std::pow(aImpact.y + bImpact.y, 2));

    vector2 n = vector2((bImpact.x - aImpact.x) / d, (bImpact.y - bImpact.y) / d);
    double p = 2 * (a.velocity * n - b.velocity * n) / 2;

    a.velocity = a.velocity - n * p * 1;
    b.velocity = b.velocity - n * p * 1;

    a.position = aImpact + a.velocity * (1 - stepProportion);
    b.position = bImpact + b.velocity * (1 - stepProportion);
}

Collision detectWallCollision(Particle p, vector2 stageSize) {
    // possible resulting particles from different choices of collision points, to pick the nearest
    vector2 end_pos = p.position + p.velocity;
    Collision result(0, 0, 2); // stepValue > 1 means no collision

    if (end_pos.x - p.radius <= 0) { // left, -1
        result = std::min(result, Collision(p.index, -1, (p.radius - p.position.x) / p.velocity.x));
    }
    if (end_pos.x + p.radius >= stageSize.x) { // right, -2
        result = std::min(result, Collision(p.index, -2, (stageSize.x - p.radius - p.position.x) / p.velocity.x));
    }
    if (end_pos.y - p.radius <= 0) { // bottom, -3
        result = std::min(result, Collision(p.index, -3, (p.radius - p.position.y) / p.velocity.y));
    }
    if (end_pos.y + p.radius >= stageSize.y) { // top, -4
        result = std::min(result, Collision(p.index, -4, (stageSize.y - p.radius - p.position.y) / p.velocity.y));
    }

    return result;
}

void resolveWallCollision(Particle& p, int wall, double stepProportion, vector2 stageSize) {
    if (wall == -1) {
        p.position.x = p.radius;
        p.position.y += stepProportion * p.velocity.y;
        p.velocity.x *= -1;
    } else if (wall == -2) {
        p.position.x = stageSize.x - p.radius;
        p.position.y += stepProportion * p.velocity.y;
        p.velocity.x *= -1;
    } else if (wall == -3) {
        p.position.x += stepProportion * p.velocity.x;
        p.position.y = p.radius;
        p.velocity.y *= -1;
    } else {
        p.position.x += stepProportion * p.velocity.x;
        p.position.y = stageSize.y - p.radius;
        p.velocity.y *= -1;
    }
}
