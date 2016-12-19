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
class LennardJonesPotential: public GlobalPotentialBase<traitsT>
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef std::pair<real_type, real_type> parameter_type;

  public:
    LennardJonesPotential() = default;
    ~LennardJonesPotential() = default;

    real_type potential(const std::size_t i, const std::size_t j,
                        const real_type r) const override;

    real_type derivative(const std::size_t i, const std::size_t j,
                         const real_type r) const override;

  private:

    std::vector<parameter_type> radii_;
};

template<typename traitsT>
inline typename LennardJonesPotential<traitsT>::real_type
LennardJonesPotential<traitsT>::potential(
        const std::size_t i, const std::size_t j, const real_type r) const
{
    const real_type r_sq   = r * r;
    const real_type sigma  = 0.5 * (radii_[i].first + radii_[j].first);
    const real_type epsilon = std::sqrt(radii_[i].second * radii_[i].second);
    const real_type sigma6_inv = 1. / std::pow(sigma, 6);
    const real_type r6s6   = r_sq * r_sq * r_sq * sigma6_inv;
    const real_type r12s12 = r6s6 * r6s6;
    return 4. * epsilon * (r12s12 - r6s6);
}

template<typename traitsT>
inline typename LennardJonesPotential<traitsT>::real_type
LennardJonesPotential<traitsT>::derivative(
        const std::size_t i, const std::size_t j, const real_type r) const
{
    const real_type r_sq   = r * r;
    const real_type r_inv  = 1. / r;
    const real_type sigma  = 0.5 * (radii_[i].first + radii_[j].first);
    const real_type epsilon = std::sqrt(radii_[i].second * radii_[i].second);
    const real_type sigma6_inv = 1. / std::pow(sigma, 6);
    const real_type r6s6   = r_sq * r_sq * r_sq * sigma6_inv;
    const real_type r12s12 = r6s6 * r6s6;

    return 24. * epsilon * r_inv * (r6s6 - 2 * r12s12);
}

} // mjolnir

#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
