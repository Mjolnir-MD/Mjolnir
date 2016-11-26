#ifndef MJOLNIR_LENNARD_JONES_POTENTIAL
#define MJOLNIR_LENNARD_JONES_POTENTIAL
#include "GlobalPotentialBase.hpp"
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
class LennardJonesPotential: public GlobalPotentialBase<traitsT,
    std::pair<typename traitsT::real_type, typename traitsT::real_type>>
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef std::pair<real_type, real_type> value_type;

  public:
    LennardJonesPotential() = default;
    ~LennardJonesPotential() = default;

    real_type
    potential(const value_type val1, const value_type val2,
              const real_type distance) const override;

    real_type
    derivative(const value_type val1, const value_type val2,
               const real_type distance) const override;
};

template<typename traitsT>
inline typename LennardJonesPotential<traitsT>::real_type
LennardJonesPotential<traitsT>::potential(
        const value_type val1, const value_type val2,
        const real_type distance) const
{
    const real_type r_sq   = distance * distance;
    const real_type sigma  = 0.5 * (val1.first + val2.first);
    const real_type epsilon = std::sqrt(val1.second + val2.second);
    const real_type sigma6_inv = 1. / std::pow(sigma, 6);
    const real_type r6s6   = r_sq * r_sq * r_sq * sigma6_inv;
    const real_type r12s12 = r6s6 * r6s6;
    return 4. * epsilon * (r12s12 - r6s6);
}

template<typename traitsT>
inline typename LennardJonesPotential<traitsT>::real_type
LennardJonesPotential<traitsT>::derivative(
        const value_type val1, const value_type val2,
        const real_type distance) const
{
    const real_type sigma  = 0.5 * (val1.first + val2.first);
    const real_type epsilon = std::sqrt(val1.second + val2.second);
    const real_type r_inv  = 1. / distance;
    const real_type r_sq   = distance;
    const real_type r6s6   = r_sq * r_sq * r_sq / std::pow(sigma, 6);
    const real_type r12s12 = r6s6 * r6s6;
    return 24. * epsilon * r_inv * (r6s6 - 2 * r12s12);
}

} // mjolnir

#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
