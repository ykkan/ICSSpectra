#ifndef ICS_RANGE_H
#define ICS_RANGE_H

#include "numeric.hpp"

struct Range{
    typedef numeric::value_type value_type;
    typedef numeric::size_type  size_type;

    value_type start;
    value_type stop;
    size_type  nstep;
    Range(): start(0), stop(0), nstep(0) {};
    Range(value_type start0, value_type stop0): start(start0), stop(stop0), nstep(1) {};
    Range(value_type start0, value_type stop0, size_type nstep0): start(start0), stop(stop0), nstep(nstep0) {};
    Range(const Range&) = default;
    ~Range() = default;
};

inline std::ostream& operator<<(std::ostream& os, const Range& r)
{
    return os << r.start << " "<< r.stop << " " << r.nstep;
}

#endif