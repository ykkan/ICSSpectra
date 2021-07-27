#include "../src/include/simulation.hpp"

int main()
{
    auto comm = Communicator();

    auto a0 = 0.5;

    auto gamma0 = 5.0;
    auto v0 = std::sqrt(1.0 - 1.0 / (gamma0 * gamma0));
    auto pp = v0 / std::sqrt(1.0 - v0 * v0);

    auto mean_x0 = Vector3{0.0, 0.0, 0.0};
    auto mean_p0 = Vector3{0.0, 0.0, pp};
    auto sigma_x0 = Vector3{1.0, 1.0, 1.0};
    auto sigma_p0 = Vector3{0.0, 0.0, 0.0};

    auto num_period = 7.0;
    auto T0 = num_period * 2 * 3.14159;
    auto beam = load_beam("beam.in");
    auto field = SinFiniteField(a0, T0, T0 / 2);

    auto dt = 0.025;
    auto nstep = static_cast<int>(T0 / dt);
    auto time_start = 0.0;
    auto time_end = time_start + dt * nstep;
    auto timewindow = Range(time_start, time_end, nstep);

    auto pi = 3.14159256;
    auto th_range = 1.0 / gamma0;
    auto th_min = pi / 2 - th_range;
    auto th_max = pi / 2 + th_range;

    auto thx_range = Range(pi / 2, pi / 2, 1);
    auto thy_range = Range(th_min, th_max, 101);
    auto t_range = Range(-2.5 * pi, 2.5 * pi, 2502);

#ifdef TRIAL
    t_range = Range(-1.070, 2.00, 3);
    auto detector =
        make_detector(thx_range, thy_range, t_range,
                      beam.num_particle, comm);
    auto tu_bound = trial_run(beam, field, detector, comm, timewindow);
#else
    auto detector = make_detector(thx_range, thy_range, t_range,
                                  beam.num_particle, comm);
    run_simulation(beam, field, detector, comm, timewindow);
    dump_detector(detector, comm, "radiation");
    comm.dump_excution_time();
#endif
}
