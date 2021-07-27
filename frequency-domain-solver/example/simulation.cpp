#include <utility>
#include <iostream>
#include "../src/include/simulation.hpp"

typedef numeric::value_type value_type;
typedef numeric::size_type size_type;

int main()
{
    Simulation sim;

    auto a0 = 0.5;
    auto num_period = 7.0;
    auto T0 = num_period * 2 * 3.141592;
    sim.set_field(field::SinFiniteEField(a0, T0, T0 / 2 + 1), field::SinFiniteBField(a0, T0, T0 / 2 + 1));

    auto beam = load_beam("beam.in");
    sim.set_beam(beam);
    sim.scatter_beam();

    auto pi = 3.14159;
    auto gamma0 = 5.0;
    auto open_angle = 1.0 / gamma0;
    auto th_min = pi / 2 - open_angle;
    auto th_max = pi / 2 + open_angle;
    auto thx_range = Range(pi / 2, pi / 2, 1);
    auto thy_range = Range(th_min, th_max, 101);
    auto w_range = Range(0.0, 500, 1251);
    sim.make_detector(thx_range, thy_range, w_range);

    auto dt = 0.025;
    auto nstep = static_cast<int>(T0 / dt);
    auto timewindow = Range(0, T0, nstep);
    sim.run(timewindow);
    sim.superimpose_radiation();
    sim.dump_detector("spectralangular.dat");

    sim.dump_excution_time();
}
