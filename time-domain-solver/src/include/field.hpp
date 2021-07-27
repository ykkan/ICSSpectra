#ifndef ICS_FIELD_H
#define ICS_FIELD_H

#include <cmath>
#include <functional>
#include "numeric.hpp"
#include "vector3.hpp"

template <typename DerivedField>
class Field
{
public:
    Vector3 efield(const Vector3 &pos, numeric::value_type t)
    {
        return static_cast<DerivedField &>(*this).efield(pos, t);
    }
    Vector3 bfield(const Vector3 &pos, numeric::value_type t)
    {
        return static_cast<DerivedField &>(*this).bfield(pos, t);
    }
};

class SinField : public Field<SinField>
{
public:
    SinField(numeric::value_type a) : amplitude(a){};
    Vector3 efield(const Vector3 &pos, numeric::value_type t)
    {
        return amplitude * Vector3(0.0, cos(t + pos.z), 0.0);
    }
    Vector3 bfield(const Vector3 &pos, numeric::value_type t)
    {
        return amplitude * Vector3(cos(t + pos.z), 0.0, 0.0);
    }

private:
    numeric::value_type amplitude;
};

numeric::value_type rect(numeric::value_type z, numeric::value_type t, numeric::value_type w, numeric::value_type z0)
{
    return abs(t + z - z0) < w / 2 ? 1.0 : 0.0;
}

class SinFiniteField : public Field<SinFiniteField>
{
public:
    SinFiniteField(numeric::value_type a, numeric::value_type width, numeric::value_type z0) : amplitude(a), width(width), z0(z0){};
    Vector3 efield(const Vector3 &pos, numeric::value_type t)
    {
        return amplitude * rect(pos.z, t, width, z0) * Vector3(0.0, cos(t + pos.z - z0), 0.0);
    }
    Vector3 bfield(const Vector3 &pos, numeric::value_type t)
    {
        return amplitude * rect(pos.z, t, width, z0) * Vector3(cos(t + pos.z - z0), 0.0, 0.0);
    }

private:
    numeric::value_type amplitude;
    numeric::value_type width;
    numeric::value_type z0;
};
#endif
