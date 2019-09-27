#pragma once

#include <cmath>

struct vector2;

double magnitude(const vector2& p);

struct vector2 {
    double x, y;

    vector2() { }

    vector2(double x, double y)
        : x(x), y(y) { };

    constexpr vector2& operator+=(const vector2& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }
    constexpr vector2& operator-=(const vector2& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }
    constexpr vector2& operator*=(const double& rhs) {
        x *= rhs;
        y *= rhs;
        return *this;
    }
    constexpr vector2& operator/=(const double& rhs) {
        x /= rhs;
        y /= rhs;
        return *this;
    }
    vector2 operator-() const {
        return { -x, -y };
    }

    void normalize()
    {
        double mag = magnitude(*this);
        x /= mag;
        y /= mag;
    }
};

inline vector2 operator+(vector2 lhs, const vector2& rhs) {
    lhs += rhs;
    return lhs;
}
inline vector2 operator-(vector2 lhs, const vector2& rhs) {
    lhs -= rhs;
    return lhs;
}
inline vector2 operator*(vector2 lhs, const double& rhs) {
    lhs *= rhs;
    return lhs;
}
inline vector2 operator/(vector2 lhs, const double& rhs) {
    lhs /= rhs;
    return lhs;
}
inline double operator*(const vector2& lhs, const vector2& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y;
}

constexpr double distSq(const vector2& p, const vector2& q) {
    return (p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y);
}
double dist(const vector2& p, const vector2& q) {
    return std::sqrt(distSq(p, q));
}
double dist(const vector2& p) {
    return dist(p, { 0, 0 });
}

inline double magnitudeSq(const vector2& p) {
    return p * p;
}

inline double magnitude(const vector2& p) {
    return std::sqrt(p * p);
}

inline double dotProduct(const vector2& a, const vector2& b)
{
    return a * b;
}
