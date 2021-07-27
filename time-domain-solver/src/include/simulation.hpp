#ifndef SIMULATION_H
#define SIMULATION_H

#include <string>
#include <random>
#include <limits>
#include <algorithm>
#include <tuple>
#include <stdlib.h>

#include "numeric.hpp"
#include "vector3.hpp"
#include "range.hpp"
#include "field.hpp"
#include "beam.hpp"
#include "detector.hpp"
#include "mpi.hpp"

typedef numeric::value_type value_type;
typedef numeric::size_type size_type;

template <typename DerivedField>
void push_beam(Beam &beam, numeric::value_type t, numeric::value_type dt, Field<DerivedField> &f)
{
    auto num_particle = beam.num_particle;
    auto q = beam.single_charge;
    auto m = beam.single_mass;
    for (Beam::size_type i = 0; i < num_particle; ++i)
    {
        auto x = beam.positions[i];
        auto p = beam.momentums[i];

        auto x_half = x + 0.5 * dt * p / sqrt(1.0 + dot(p, p));
        auto p_minus = p + 0.5 * (q / m) * dt * f.efield(x_half, t + 0.5 * dt);

        auto gamma_minus = sqrt(1.0 + dot(p_minus, p_minus));
        auto t_vec = 0.5 * dt * (q / m) * f.bfield(x_half, t + 0.5 * dt) / gamma_minus;
        auto s_vec = 2.0 * t_vec / (1.0 + dot(t_vec, t_vec));

        auto p_plus = p_minus + cross(p_minus + (cross(p_minus, t_vec)), s_vec);

        auto p_final = p_plus + 0.5 * dt * (q / m) * f.efield(x_half, t + 0.5 * dt);
        auto v_final = p_final / sqrt(1.0 + dot(p_final, p_final));

        auto x_final = x_half + 0.5 * dt * v_final;

        auto gamma_final = sqrt(1.0 + dot(p_final, p_final));

        beam.positions[i] = x_final;
        beam.momentums[i] = p_final;
        beam.velocities[i] = v_final;

        auto F = q * (f.efield(x_final, t + dt) + cross(v_final, f.bfield(x_final, t + dt)));
        beam.accelerations[i] = 1 / (m * gamma_final) * (F - dot(v_final, F) * v_final);
    }
}

void recieve_radiation(Detector &detector, Beam &beam, value_type t)
{
    for (auto &s : detector)
    {
        s.recieve_field(beam, t);
    }
}


Detector make_detector(const Range &thx_range,
                       const Range &thy_range,
                       const Range &t_range,
                       size_type npar, Communicator &communicator)
{
    auto tx_min = thx_range.start;
    auto tx_max = thx_range.stop;
    auto ntx = thx_range.nstep;

    auto ty_min = thy_range.start;
    auto ty_max = thy_range.stop;
    auto nty = thy_range.nstep;

    auto t_begin = t_range.start;
    auto t_end = t_range.stop;
    auto ntsample = t_range.nstep;

    ntx = (ntx % 2 == 0) ? (ntx + 1) : ntx;
    nty = (nty % 2 == 0) ? (nty + 1) : nty;
    auto dtx = (tx_max - tx_min) / (ntx - 1);
    auto dty = (ty_max - ty_min) / (nty - 1);

    auto nprocs = communicator.size();
    auto procs_id = communicator.rank();

    auto nsensor = ntx * nty;
    auto nsensor_local = (procs_id < (nsensor % nprocs)) ? (nsensor / nprocs + 1) : (nsensor / nprocs);

    auto sensor_default = Sensor(Vector3(0, 0, 0), npar, t_begin, t_end, ntsample);
    std::vector<Sensor> sensors(nsensor_local, sensor_default);

    auto start_index = (procs_id <= nsensor % nprocs) ? (nsensor / nprocs + 1) * procs_id : (nsensor / nprocs + 1) * (nsensor % nprocs) + (nsensor / nprocs) * (procs_id - nsensor % nprocs);

    for (auto i = 0; i < nsensor_local; ++i)
    {
        auto sensor_index = start_index + 1 + i;
        auto tx_index = sensor_index / nty;
        auto ty_index = sensor_index % nty;
        auto tx = (ntx > 1) ? (tx_min + tx_index * dtx) : tx_min;
        auto ty = (nty > 1) ? (ty_min + ty_index * dty) : ty_min;
        sensors[i].position = Vector3(cos(tx), sin(tx) * cos(ty), sin(tx) * sin(ty));
        sensors[i].index = sensor_index;
    }

    //dump detector info
    if (procs_id == 0)
    {
        std::ofstream file;
        auto filename = "detectorinfo.dat";
        file.open(filename);
        file << "thetaXRange: " << tx_min << " " << tx_max << " " << ntx << "\n";
        file << "thetaYrange: " << ty_min << " " << ty_max << " " << nty << "\n";
        file << "pulseRange: " << t_begin << " " << t_end << " " << ntsample << "\n";
        file.close();
    }
    return sensors;
}

void dump_detector(const Detector &detector, Communicator &communicator, const std::string &file_prefix)
{
    if (communicator.rank() == 0)
    {
        //system("rm -rf " + file_prefix);
        system("mkdir -p rawdata");
    }
    communicator.barrier();

    for (auto &s : detector)
    {
        std::ofstream file;
        auto filename = file_prefix + "_" + std::to_string(s.index) + ".dat";
        file.open("rawdata/" + filename);
        for (auto &f : s.total_fields)
        {
            file << f << "\n";
        }
        file.close();
    }
    communicator.barrier();
}

auto beam_ta_min = [](Beam &beam, Vector3 &n_vec, value_type t)
{
    auto ta_min = std::numeric_limits<value_type>::max();
    for (auto &x : beam.positions)
    {
        auto ta = t - dot(n_vec, x);
        ta_min = (ta_min < ta) ? ta_min : ta;
    }
    return ta_min;
};

auto beam_ta_max = [](Beam &beam, Vector3 &n_vec, value_type t)
{
    auto ta_max = std::numeric_limits<value_type>::min();
    for (auto &x : beam.positions)
    {
        auto ta = t - dot(n_vec, x);
        ta_max = (ta_max > ta) ? ta_max : ta;
    }
    return ta_max;
};


template <typename DerivedField>
Range trial_run(Beam &beam, Field<DerivedField> &field, Detector &detector, Communicator &comm,
                const Range &time_window)
{
    auto t_begin = time_window.start;
    auto t_end = time_window.stop;
    auto nstep = time_window.nstep;
    auto dt = (t_end - t_begin) / nstep;

    auto ta_begin_min = std::numeric_limits<value_type>::max();
    for (auto &s : detector)
    {
        auto ta_min = beam_ta_min(beam, s.position, t_begin);
        ta_begin_min = (ta_begin_min < ta_min) ? ta_begin_min : ta_min;
    }
    for (int i = 1; i <= nstep; ++i)
    {
        if (comm.rank() == 0)
        {
            std::cout << "trial step: " << i << "/ " << nstep << "\n";
        }
        auto t = t_begin + i * dt;
        push_beam(beam, t, dt, field);
    }
    auto ta_end_max = std::numeric_limits<value_type>::min();
    for (auto &s : detector)
    {
        auto ta_max = beam_ta_max(beam, s.position, t_end);
        ta_end_max = ta_end_max > ta_max ? ta_end_max : ta_max;
    }

    auto ta_min = ta_begin_min;
    auto ta_max = ta_end_max;
    value_type ta_min_global = 0;
    value_type ta_max_global = 0;
    comm.find_min(ta_min, ta_min_global);
    comm.find_max(ta_max, ta_max_global);

    if (comm.rank() == 0)
    {
        std::ofstream file;
        auto filename = "trialrun.dat";
        file.open(filename);
        file << "result of trial run\n";
        file << "ta begin: " << ta_min_global << "\n";
        file << "ta end  : " << ta_max_global << "\n";
        file << "\n";
        file.close();

        std::cout << "result of trial run\n";
        std::cout << "ta begin: " << ta_min_global << "\n";
        std::cout << "ta end  : " << ta_max_global << "\n";
        std::cout << "\n";
    }
    return Range(ta_min_global, ta_max_global);
    comm.barrier();
}

template <typename DerivedField>
void run_simulation(Beam &beam, Field<DerivedField> &field, Detector &detector, Communicator &comm,
                    const Range &time_window)
{
    auto t_begin = time_window.start;
    auto t_end = time_window.stop;
    auto nstep = time_window.nstep;
    auto dt = (t_end - t_begin) / nstep;

    for (int i = 0; i < nstep; ++i)
    {
        if (comm.rank() == 0)
        {
            std::cout << "run step: " << i << "/ " << nstep << "\n";
        }
        auto t = t_begin + i * dt;
        recieve_radiation(detector, beam, t);
        push_beam(beam, t, dt, field);
    }
}

#endif