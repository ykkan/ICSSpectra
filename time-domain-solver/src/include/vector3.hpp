#ifndef ICS_VECTOR3_H
#define ICS_VECTOR3_H

#include <iostream>
#include "numeric.hpp"

struct Vector3
{
    typedef numeric::value_type value_type;

    value_type x;
    value_type y;
    value_type z;

    Vector3() : x(0), y(0), z(0){};
    Vector3(value_type x0, value_type y0, value_type z0) : x(x0), y(y0), z(z0){};
    Vector3(const Vector3 &) = default;
    ~Vector3() = default;
    Vector3 &operator=(const Vector3 &) = default;
    Vector3 &operator+=(const Vector3 &v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
        return *this;
    }
    Vector3 &operator-=(const Vector3 &v)
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
        return *this;
    }
    Vector3 &operator*=(const value_type &c)
    {
        x *= c;
        y *= c;
        z *= c;
        return *this;
    }
    Vector3 &operator/=(const value_type &c)
    {
        x /= c;
        y /= c;
        z /= c;
        return *this;
    }
    bool operator==(const Vector3 &other) const
    {
        return (this->x == other.x && this->y == other.y && this->z == other.z);
    }
    bool operator!=(const Vector3 &other) const
    {
        return !(*this == other);
    }
};

inline Vector3 operator*(const Vector3::value_type &c, const Vector3 &v)
{
    return (Vector3(v) *= c);
}

inline Vector3 operator*(const Vector3 &v, const Vector3::value_type &c)
{
    return (Vector3(v) *= c);
}

inline Vector3 operator/(const Vector3 &v, const Vector3::value_type &c)
{
    return (Vector3(v) /= c);
}

inline Vector3 operator+(const Vector3 &v1, const Vector3 &v2)
{
    return (Vector3(v1) += v2);
}

inline Vector3 operator-(const Vector3 &v1, const Vector3 &v2)
{
    return (Vector3(v1) -= v2);
}

inline Vector3::value_type dot(const Vector3 &v1, const Vector3 &v2)
{
    return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}

inline Vector3 cross(const Vector3 &v1, const Vector3 &v2)
{
    return Vector3(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

inline Vector3::value_type norm2(const Vector3 &v)
{
    return (v.x * v.x + v.y * v.y + v.z * v.z);
}

inline std::ostream &operator<<(std::ostream &os, const Vector3 &v)
{
    return os << v.x << " " << v.y << " " << v.z;
}

#endif