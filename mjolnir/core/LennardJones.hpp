#ifndef MJOLNIR_LENNARD_JONES_POTENTIAL
#define MJOLNIR_LENNARD_JONES_POTENTIAL
#include <mjolnir/math/Vector.hpp>
#include <mjolnir/math/fast_inv_sqrt.hpp>
#include <cmath>

namespace mjolnir
{

/*! @brief Lennard-Jones type potential & derivative                       *
 * designed for global force field. so it doesn't have its own parameters. *
 * V(r)  =  4. * epsilon * ((r/sigma)^12 - (r/sigma)^6))                   *
 * dV/dr = 24. * epsilon / r * ((r/sigma)^6 - 2 * (r/sigma)^12)            */
template<typename traitsT>
class LennardJones
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:
    LennardJones() = default;
    ~LennardJones() = default;

    real_type
    potential(const real_type length_square,
              const real_type sigma, const real_type epsilon) const;

    real_type
    derivative(const real_type length_square,
               const real_type sigma, const real_type epsilon) const;
};

template<typename traitsT>
inline typename LennardJones<traitsT>::real_type
LennardJones<traitsT>::potential(const real_type length_square,
              const real_type sigma, const real_type epsilon) const
{
    const real_type r_sq   = length_square;
    const real_type sigma6_inv = 1. / std::pow(sigma, 6);
    const real_type r6s6   = r_sq * r_sq * r_sq * sigma6_inv;
    const real_type r12s12 = r6s6 * r6s6;
    return 4. * epsilon * (r12s12 - r6s6);
}

template<typename traitsT>
inline typename LennardJones<traitsT>::real_type
LennardJones<traitsT>::derivative(const real_type length_square,
              const real_type sigma, const real_type epsilon) const
{
    const real_type r_inv  = fast_inv_sqrt(length_square);
    const real_type r_sq   = length_square;
    const real_type r6s6   = r_sq * r_sq * r_sq / std::pow(sigma, 6);
    const real_type r12s12 = r6s6 * r6s6;
    return 24. * epsilon * r_inv * (r6s6 - 2 * r12s12);
}

} // mjolnir

#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
