#ifndef ICS_BEAM_H
#define ICS_BEAM_H

#include <cmath>
#include <vector>
#include <string>
#include <fstream>

#include "numeric.hpp"
#include "vector3.hpp"
#include "field.hpp"

struct Beam
{
    typedef numeric::value_type value_type;
    typedef numeric::size_type size_type;
    typedef std::vector<Vector3> container_type;

    size_type num_particle;
    value_type single_charge;
    value_type single_mass;
    container_type positions;
    container_type momentums;

    Beam() = default;

    Beam(size_type npar) : num_particle(npar),
                           single_charge(-1),
                           single_mass(1),
                           positions(npar), momentums(npar){};

    Beam(size_type npar, value_type c, value_type m) : num_particle(npar),
                                                       single_charge(c),
                                                       single_mass(m),
                                                       positions(npar), momentums(npar){};

    Beam(size_type npar, value_type c, value_type m, const Vector3 &x0, const Vector3 &p0) : num_particle(npar),
                                                                                             single_charge(c),
                                                                                             single_mass(m),
                                                                                             positions(npar, x0), momentums(npar, p0){};
    Beam(Beam &&) = default;

    Beam &operator=(Beam &&) = default;
};

Beam load_beam(const std::string &filename)
{
    std::ifstream infile(filename);
    numeric::value_type x;
    numeric::value_type y;
    numeric::value_type z;
    numeric::value_type px;
    numeric::value_type py;
    numeric::value_type pz;
    numeric::value_type charge;
    numeric::value_type mass;
    numeric::size_type npar;

    infile >> npar >> charge >> mass;
    auto beam = Beam(npar, charge, mass);
    for (numeric::size_type i = 0; i < npar; ++i)
    {
        infile >> x >> y >> z >> px >> py >> pz;
        beam.positions[i] = Vector3{x, y, z};
        beam.momentums[i] = Vector3{px, py, pz};
    }
    return beam;
}

#endif
