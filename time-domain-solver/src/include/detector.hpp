#ifndef ICS_DETECTOR_H
#define ICS_DETECTOR_H

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <tuple>

#include "numeric.hpp"
#include "vector3.hpp"
#include "range.hpp"
#include "field.hpp"

struct Sensor
{
    typedef numeric::value_type value_type;
    typedef numeric::size_type size_type;

    Vector3 position;

    value_type num_particle;

    value_type t_begin;
    value_type t_end;
    value_type dt;
    size_type num_tsample;

    std::vector<value_type> buff_ts;
    std::vector<Vector3> buff_fields;

    bool is_interpolatable = false;
    std::vector<Vector3> total_fields;

    size_type index;

    Sensor(Vector3 pos, size_type npar, value_type t_begin, value_type t_end, size_type num_tsample) : position(pos),
                                                                                                       num_particle(npar),
                                                                                                       t_begin(t_begin),
                                                                                                       t_end(t_end),
                                                                                                       num_tsample(num_tsample),
                                                                                                       dt((t_end - t_begin) / (num_tsample - 1)),
                                                                                                       buff_ts(npar),
                                                                                                       buff_fields(npar),
                                                                                                       total_fields(num_tsample),
                                                                                                       index(0)
    {
    }
    void recieve_field(const Beam &beam, value_type t)
    {
        auto n_vec = position;
        auto q = beam.single_charge;
        if (is_interpolatable)
        {
            for (size_type k = 0; k < num_particle; ++k)
            {
                auto v = beam.velocities[k];
                auto v_dot = beam.accelerations[k];

                auto ta_pre = buff_ts[k];
                auto ta_now = t - dot(n_vec, beam.positions[k]);

                auto field_pre = buff_fields[k];
                auto field_now = q * cross(n_vec, cross(n_vec - v, v_dot)) / std::pow((1.0 - dot(n_vec, v)), 3);

                // determine the upper and lower indecies of interpolation points
                // ceil
                auto index_begin = static_cast<size_type>((ta_pre - t_begin) / dt) + 1; 
                // floor
                auto index_end = static_cast<size_type>((ta_now - t_begin) / dt);
                      
                for (auto i = index_begin; i <= index_end; ++i)
                {
                    // linear interpolation
                    auto ti = t_begin + dt * i;
                    total_fields[i] += (field_pre * (ta_now - ti) + field_now * (ti - ta_pre)) / (ta_now - ta_pre);
                }

                buff_ts[k] = ta_now;
                buff_fields[k] = field_now;
            }
        }
        else
        {
            for (size_type i = 0; i < num_particle; ++i)
            {
                auto &v = beam.velocities[i];
                auto &v_dot = beam.accelerations[i];
                buff_ts[i] = t - dot(n_vec, beam.positions[i]);
                buff_fields[i] = q * cross(n_vec, cross(n_vec - v, v_dot)) / std::pow((1.0 - dot(n_vec, v)), 3);
            }
            is_interpolatable = true;
        }
    }
};

using Detector = std::vector<Sensor>;

#endif
