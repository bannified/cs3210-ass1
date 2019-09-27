#pragma once

#include "vector2.h"
#include "Particle.h"

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
