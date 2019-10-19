#pragma once

#include <cmath>

struct vector2;

__device__ double magnitude(const vector2& p);

struct vector2
{
    double x, y;

    vector2()
    {
    }

    vector2(double x, double y)
        : x(x), y(y)
    {
    };

    __device__ vector2& operator+=(const vector2& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }
    __device__ vector2& operator-=(const vector2& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }
    __device__ vector2& operator*=(const double& rhs)
    {
        x *= rhs;
        y *= rhs;
        return *this;
    }
    __device__ vector2& operator/=(const double& rhs)
    {
        x /= rhs;
        y /= rhs;
        return *this;
    }
    __device__ vector2 operator-() const
    {
        return { -x, -y };
    }

    __device__ void normalize()
    {
        double mag = magnitude(*this);
        x /= mag;
        y /= mag;
    }
};

__device__ inline vector2 operator+(vector2 lhs, const vector2& rhs)
{
    lhs += rhs;
    return lhs;
}
__device__ inline vector2 operator-(vector2 lhs, const vector2& rhs)
{
    lhs -= rhs;
    return lhs;
}
__device__ inline vector2 operator*(vector2 lhs, const double& rhs)
{
    lhs *= rhs;
    return lhs;
}
__device__ inline vector2 operator/(vector2 lhs, const double& rhs)
{
    lhs /= rhs;
    return lhs;
}
__device__ inline double operator*(const vector2& lhs, const vector2& rhs)
{
    return lhs.x * rhs.x + lhs.y * rhs.y;
}

__device__ double distSq(const vector2& p, const vector2& q)
{
    return (p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y);
}
__device__ double dist(const vector2& p, const vector2& q)
{
    return std::sqrt(distSq(p, q));
}

__device__ double dist(const vector2& p)
{
    return std::sqrt((p.x) * (p.x) + (p.y) * (p.y));
}

__device__ double magnitudeSq(const vector2& p)
{
    return p * p;
}

__device__ double magnitude(const vector2& p)
{
    return std::sqrt(p * p);
}

__device__ double dotProduct(const vector2& a, const vector2& b)
{
    return a * b;
}
