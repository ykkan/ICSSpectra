#ifndef ICS_FIELD_H
#define ICS_FIELD_H

#include <cmath>
#include <functional>
#include "numeric.hpp"
#include "vector3.hpp"

namespace field
{

    typedef numeric::value_type value_type;
    typedef std::function<Vector3(Vector3 x, value_type t)> Field;

    class SinEField
    {
    public:
        SinEField(value_type a) : amplitude(a){};
        Vector3 operator()(Vector3 pos, numeric::value_type t) const
        {
            return amplitude * Vector3(0.0, cos(t + pos.z), 0.0);
        }

    private:
        numeric::value_type amplitude;
    };

    class SinBField
    {
    public:
        SinBField(value_type a) : amplitude(a){};
        Vector3 operator()(Vector3 pos, numeric::value_type t) const
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

    class SinFiniteEField
    {
    public:
        SinFiniteEField(value_type a, value_type w, value_type z0) : amplitude(a), w(w), z0(z0){};
        Vector3 operator()(Vector3 pos, numeric::value_type t) const
        {
            return amplitude * rect(pos.z, t, w, z0) * Vector3(0.0, cos(t + pos.z - z0), 0.0);
        }

    private:
        numeric::value_type amplitude;
        numeric::value_type w;
        numeric::value_type z0;
    };

    class SinFiniteBField
    {
    public:
        SinFiniteBField(value_type a, value_type w, value_type z0) : amplitude(a), w(w), z0(z0){};
        Vector3 operator()(Vector3 pos, numeric::value_type t) const
        {
            return amplitude * rect(pos.z, t, w, z0) * Vector3(cos(t + pos.z - z0), 0.0, 0.0);
        }

    private:
        numeric::value_type amplitude;
        numeric::value_type w;
        numeric::value_type z0;
    };

    class GaussianEField
    {
    public:
        GaussianEField(value_type a, value_type waist, value_type duration) : amplitude(a),
                                                                              W0(waist),
                                                                              tau(duration),
                                                                              ZR(0.5 * waist * waist){};
        Vector3 operator()(Vector3 pos, numeric::value_type t) const
        {
            auto x = pos.x;
            auto y = pos.y;
            auto z = pos.z;
            auto r = sqrt(x * x + y * y);

            auto R = z + ZR * ZR / z;
            auto W = W0 * sqrt(1.0 + pow(z / ZR, 2.0));

            auto Ex = amplitude / sqrt(1.0 + pow(z / ZR, 2.0)) *
                      exp(-pow(r / W, 2.0) - 1.38629436112 * pow(t / tau, 2.0)) *
                      cos(t + z - r * r / (2 * R) + atan(z / ZR));
            return Vector3(amplitude * Ex, 0.0, 0.0);
        }

    private:
        value_type amplitude;
        value_type W0;
        value_type tau;
        value_type ZR;
    };

    class GaussianBField
    {
    public:
        GaussianBField(Field ef)
        {
            efield = ef;
        }
        Vector3 operator()(Vector3 pos, numeric::value_type t) const
        {
            auto ef = efield(pos, t);
            return Vector3(0.0, -ef.x, 0.0);
        }

    private:
        Field efield;
    };

}

/*
Vector3 efield(const Vector3& pos, numeric::value_type t)
{
    return Vector3(
        0.0,
        cos(t + pos.z),
        0.0
    );
}

Vector3 bfield(const Vector3& pos, numeric::value_type t)
{
    return Vector3(
        cos(t + pos.z),
        0.0,
        0.0
    );
}
*/

#endif
