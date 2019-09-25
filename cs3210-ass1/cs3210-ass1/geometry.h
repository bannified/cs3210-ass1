#include "point.h"

// find closest point on line segment pq to point r
point closestPointOnLine(point p, point q, point r){
    auto [x1, y1] = p;
    auto [x2, y2] = q;
    auto [x0, y0] = r;
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
