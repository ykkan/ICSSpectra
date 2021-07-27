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
#include "beam.hpp"
#include "field.hpp"

struct Detector
{
    /* data */
    typedef numeric::value_type value_type;
    typedef numeric::size_type size_type;
    typedef std::vector<Vector3> container_type;

    value_type tx_min;
    value_type tx_max;
    size_type num_tx;

    value_type ty_min;
    value_type ty_max;
    size_type num_ty;

    value_type w_min;
    value_type w_max;
    size_type num_omega;

    value_type lambda0;

    std::vector<value_type> txs;
    std::vector<value_type> tys;
    std::vector<value_type> omegas;

    field::Field efield;
    field::Field bfield;

    container_type re_data;
    container_type im_data;
    Detector() = default;

    Detector(const Range &thx_range,
             const Range &thy_range,
             const Range &w_range,
             field::Field ef, field::Field bf) : tx_min(thx_range.start),
                                                 tx_max(thx_range.stop),
                                                 num_tx((thx_range.nstep % 2 == 0) ? thx_range.nstep + 1 : thx_range.nstep),
                                                 ty_min(thy_range.start),
                                                 ty_max(thy_range.stop),
                                                 num_ty((thy_range.nstep % 2 == 0) ? thy_range.nstep + 1 : thy_range.nstep),
                                                 w_min(w_range.start),
                                                 w_max(w_range.stop),
                                                 num_omega(w_range.nstep),
                                                 efield(ef),
                                                 bfield(bf)
    {
        txs = std::vector<value_type>(num_tx);
        tys = std::vector<value_type>(num_ty);
        omegas = std::vector<value_type>(num_omega);
        re_data = container_type(num_tx * num_ty * num_omega);
        im_data = container_type(num_tx * num_ty * num_omega);

        for (size_type i = 0; i < num_tx; ++i)
        {
            if (num_tx == 1)
            {
                txs[i] = tx_min;
            }
            else
            {
                txs[i] = tx_min + i * (tx_max - tx_min) / (num_tx - 1);
            }
        }
        for (size_type i = 0; i < num_ty; ++i)
        {
            if (num_ty == 1)
            {
                tys[i] = ty_min;
            }
            else
            {
                tys[i] = ty_min + i * (ty_max - ty_min) / (num_ty - 1);
            }
        }
        for (size_type i = 0; i < num_omega; ++i)
        {
            if (num_omega == 1)
            {
                omegas[i] = w_min;
            }
            else
            {
                omegas[i] = w_min + i * (w_max - w_min) / (num_omega - 1);
            }
        }
    }

    void accumulate_radiation(const Beam &beam, numeric::value_type t, numeric::value_type dt);

    void dump_data(const std::string &filename);

private:
    // map sensor 2d coordinate to 1d index
    size_type coord2index(const size_type &i, const size_type &j)
    {
        return (j + i * num_tx);
    }
    void update_bunch_field(const Beam &beam, const size_type &sensor_index, const Vector3 &n_vec,
                            const value_type &w, const numeric::value_type &t, const numeric::value_type &dt);
};

using std::cos;
using std::sin;
using std::sqrt;

void Detector::update_bunch_field(const Beam &beam, const size_type &sensor_index,
                                  const Vector3 &n_vec, const value_type &w,
                                  const numeric::value_type &t, const numeric::value_type &dt)
{
    auto num_particle = beam.num_particle;
    auto q = beam.single_charge;
    auto m = beam.single_mass;

    auto &re_field = re_data[sensor_index];
    auto &im_field = im_data[sensor_index];

    for (numeric::size_type i = 0; i < num_particle; ++i)
    {
        auto x = beam.positions[i];
        auto p = beam.momentums[i];

        auto gamma = sqrt(1.0 + dot(p, p));
        auto v = p / gamma;

        auto F = q * (efield(x, t) + cross(v, bfield(x, t)));
        auto v_dot = 1.0 / (m * gamma) * (F - dot(v, F) * v);
        auto temp_field = q * dt * cross(n_vec, cross(n_vec - v, v_dot)) / std::pow((1.0 - dot(n_vec, v)), 2);
        auto phase = w * (t - dot(n_vec, x));

        re_field += temp_field * cos(phase);
        im_field += temp_field * sin(phase);
    }
}

void Detector::accumulate_radiation(const Beam &beam, numeric::value_type t, numeric::value_type dt)
{
    for (numeric::size_type i = 0; i < num_tx; ++i)
    {
        for (numeric::size_type j = 0; j < num_ty; ++j)
        {
            auto tx = txs[i];
            auto ty = tys[j];
            auto n_vec = Vector3(cos(tx), sin(tx) * cos(ty), sin(tx) * sin(ty));
            auto sensor_start_index = num_omega * coord2index(i, j);
            for (numeric::size_type k = 0; k < num_omega; ++k)
            {
                update_bunch_field(beam, sensor_start_index + k, n_vec, omegas[k], t, dt);
            }
        }
    }
}

void Detector::dump_data(const std::string &filename)
{

    std::ofstream file;
    file.open(filename);
    for (numeric::size_type i = 0; i < num_tx; ++i)
    {
        for (numeric::size_type j = 0; j < num_ty; ++j)
        {
            auto sensor_start_index = num_omega * coord2index(i, j);
            for (numeric::size_type k = 0; k < num_omega; ++k)
            {
                file << norm2(re_data[sensor_start_index + k]) + norm2(im_data[sensor_start_index + k]) << " ";
            }
            file << "\n";
        }
    }
    file.close();
}
#endif