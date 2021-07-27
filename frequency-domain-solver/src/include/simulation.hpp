#ifndef SIMULATION_H
#define SIMULATION_H

#include <string>
#include <random>

#include "numeric.hpp"
#include "vector3.hpp"
#include "range.hpp"
#include "field.hpp"
#include "beam.hpp"
#include "detector.hpp"
#include "mpi.hpp"

class Simulation
{
public:
    typedef numeric::value_type value_type;
    typedef numeric::size_type size_type;

    //Simulation() = default;
    void set_beam(Beam &beam);
    void scatter_beam();

    void push_beam(numeric::value_type t, numeric::value_type dt);

    void set_field(field::Field efield, field::Field bfield);

    void run(const Range &time_window);

    void make_detector(const Range &thx_range,
                       const Range &thy_range,
                       const Range &w_range);

    void superimpose_radiation();

    void dump_detector(const std::string &filename);

    void dump_excution_time();

public:
    Communicator comm;

    Beam::size_type num_particle_global;
    Beam::size_type num_particle_local;
    Beam beam_global;
    Beam beam_local;

    field::Field efield;
    field::Field bfield;

    Detector detector_global;
    Detector detector_local;
};

void Simulation::set_field(field::Field ef, field::Field bf)
{
    efield = ef;
    bfield = bf;
}

void Simulation::push_beam(numeric::value_type t, numeric::value_type dt)
{
    auto num_particle = beam_local.num_particle;
    auto q = beam_local.single_charge;
    auto m = beam_local.single_mass;
    for (Beam::size_type i = 0; i < num_particle; ++i)
    {
        auto x = beam_local.positions[i];
        auto p = beam_local.momentums[i];

        auto x_half = x + 0.5 * dt * p / sqrt(1.0 + dot(p, p));
        auto p_minus = p + 0.5 * (q / m) * dt * efield(x_half, t + 0.5 * dt);

        auto gamma_minus = sqrt(1.0 + dot(p_minus, p_minus));
        auto t_vec = 0.5 * dt * (q / m) * bfield(x_half, t + 0.5 * dt) / gamma_minus;
        auto s_vec = 2.0 * t_vec / (1.0 + dot(t_vec, t_vec));

        auto p_plus = p_minus + cross(p_minus + (cross(p_minus, t_vec)), s_vec);

        auto p_final = p_plus + 0.5 * dt * (q / m) * efield(x_half, t + 0.5 * dt);
        auto v_final = p_final / sqrt(1.0 + dot(p_final, p_final));

        auto x_final = x_half + 0.5 * dt * v_final;

        beam_local.positions[i] = x_final;
        beam_local.momentums[i] = p_final;
    }
}

void Simulation::set_beam(Beam &beam)
{
    num_particle_global = beam.num_particle;
    if (comm.rank() == 0)
    {
        beam_global = std::move(beam);
    }
}

void Simulation::scatter_beam()
{
    num_particle_local = num_particle_global / comm.size();
    num_particle_local += comm.rank() < (num_particle_global % comm.size()) ? 1 : 0;

    //initalize information for local beam
    beam_local = Beam(num_particle_local);
    value_type m;
    value_type q;
    if (comm.rank() == 0)
    {
        m = beam_global.single_mass;
        q = beam_global.single_charge;
    }
    comm.broadcast(m);
    comm.broadcast(q);

    beam_local.single_mass = m;
    beam_local.single_charge = q;

    // scatter position and momentum to local beam
    std::vector<size_type> counts(comm.size());
    for (int i = 0; i < comm.size(); ++i)
    {
        counts[i] = num_particle_global / comm.size() + (i < num_particle_global % comm.size() ? 1 : 0);
    }
    std::vector<size_type> displs(comm.size());
    int sum = 0;
    for (auto i = 0; i < comm.size(); ++i)
    {
        displs[i] = sum;
        sum += counts[i];
    }

    comm.scatter_Vector3s(beam_global.positions, beam_local.positions, counts, displs);
    comm.scatter_Vector3s(beam_global.momentums, beam_local.momentums, counts, displs);
}

void Simulation::run(const Range &time_window)
{
    auto t_begin = time_window.start;
    auto t_end = time_window.stop;
    auto nstep = time_window.nstep;
    auto dt = (t_end - t_begin) / nstep;
    for (size_type i = 0; i < nstep; ++i)
    {
        if (comm.rank() == 0)
        {
            std::cout << "step: " << i << "/ " << nstep << "\n";
        }
        auto t = t_begin + i * dt;
        detector_local.accumulate_radiation(beam_local, t, dt);
        push_beam(t, dt);

        //std::cout << t << " " << beam_local.positions[0] << "\n";
    }
}

void Simulation::make_detector(const Range &thx_range,
                               const Range &thy_range,
                               const Range &w_range)
{
    detector_global = Detector(thx_range, thy_range, w_range, efield, bfield);
    detector_local = Detector(thx_range, thy_range, w_range, efield, bfield);
}

void Simulation::superimpose_radiation()
{
    comm.sum_Vector3s(detector_local.re_data, detector_global.re_data);
    comm.sum_Vector3s(detector_local.im_data, detector_global.im_data);
}

void Simulation::dump_detector(const std::string &filename)
{
    if (comm.rank() == 0)
    {
        detector_global.dump_data(filename);
    }
}

void Simulation::dump_excution_time()
{
    comm.dump_excution_time();
}

#endif