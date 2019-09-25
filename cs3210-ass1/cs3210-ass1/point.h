#include <cmath>

struct point {
    double x, y;

    constexpr point& operator+=(const point& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }
    constexpr point& operator-=(const point& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }
    constexpr point& operator*=(const double& rhs) {
        x *= rhs;
        y *= rhs;
        return *this;
    }
    constexpr point& operator/=(const double& rhs) {
        x /= rhs;
        y /= rhs;
        return *this;
    }
    constexpr point operator-() const {
        return { -x, -y };
    }
};

constexpr inline point operator+(point lhs, const point& rhs) {
    lhs += rhs;
    return lhs;
}
constexpr inline point operator-(point lhs, const point& rhs) {
    lhs -= rhs;
    return lhs;
}
constexpr inline point operator*(point lhs, const double& rhs) {
    lhs *= rhs;
    return lhs;
}
constexpr inline point operator/(point lhs, const double& rhs) {
    lhs /= rhs;
    return lhs;
}
constexpr inline double operator*(const& point lhs, const point& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y;
}

constexpr double distSq(const point& p, const point& q) {
    return (p.x - q.x) * (p.x - q.x) + (p.y - q.y) * (p.y - q.y);
}
double dist(const point& p, const point& q) {
    return std::sqrt(distSq(p, q));
}
double dist(const point& p) {
    return dist(p, { 0, 0 });
}
