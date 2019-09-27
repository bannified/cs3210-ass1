#pragma once

#include "point.h"

// find closest point on line segment pq to point r
point closestPointOnLine(point p, point q, point r){
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
