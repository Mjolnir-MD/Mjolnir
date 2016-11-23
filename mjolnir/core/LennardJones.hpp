#ifndef MJOLNIR_LENNARD_JONES_POTENTIAL
#define MJOLNIR_LENNARD_JONES_POTENTIAL
#include <mjolnir/math/Vector.hpp>
#include <mjolnir/math/fast_inv_sqrt.hpp>
#include "PotentialBase.hpp"
#include <cmath>

namespace mjolnir
{

template<typename traitsT>
class LennardJones : public PotentialBase<traitsT, 2>
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:
    LennardJones(const real_type sigma, const real_type epsilon)
        : sigma6_inv_(1. / std::pow(sigma, 6)), epsilon_(epsilon)
    {}

    real_type
    potential(std::array<coordinate_type const&, 2> pos) const override;
    real_type
    derivative(std::array<coordinate_type const&, 2> pos) const override;

  private:

    real_type sigma6_inv_;
    real_type epsilon_;
};

template<typename traitsT>
inline typename LennardJones<traitsT>::real_type
LennardJones<traitsT>::potential(
        std::array<coordinate_type const&, 2> pos) const
{
    const real_type r_sq   = length_sq(pos[0] - pos[1]);
    const real_type r6s6   = r_sq * r_sq * r_sq * sigma6_inv_;
    const real_type r12s12 = r6s6 * r6s6;
    return 4. * epsilon_ * (r12s12 - r6s6);
}

template<typename traitsT>
inline typename LennardJones<traitsT>::real_type
LennardJones<traitsT>::derivative(
        std::array<coordinate_type const&, 2> pos) const
{
    const real_type r_inv  = fast_inv_sqrt(length_sq(pos[0] - pos[1]));
    const real_type r_sq   = length_sq(pos[0] - pos[1]);
    const real_type r6s6   = r_sq * r_sq * r_sq * sigma6_inv_;
    const real_type r12s12 = r6s6 * r6s6;
    return 24. * epsilon_ * r_inv * (r6s6 - 2 * r12s12);
}

} // mjolnir

#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
