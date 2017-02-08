#ifndef MJOLNIR_BOND_ANGLE_INTERACTION
#define MJOLNIR_BOND_ANGLE_INTERACTION
#include "LocalInteractionBase.hpp"
#include "BoundaryCondition.hpp"
#include "constants.hpp"
#include <mjolnir/math/fast_inv_sqrt.hpp>
#include <cmath>

namespace mjolnir
{

template<typename traitsT, typename boundaryT = UnlimitedBoundary<traitsT>>
class BondAngleInteraction : public LocalInteractionBase<traitsT, 3>
{
  public:
    typedef traitsT   traits_type;
    typedef boundaryT boundary_type;
    typedef LocalInteractionBase<traits_type, 3>    base_type;
    typedef typename base_type::time_type           time_type;
    typedef typename base_type::real_type           real_type;
    typedef typename base_type::coordinate_type     coordinate_type;
    typedef typename base_type::particle_type       particle_type;
    typedef LocalPotentialBase<traits_type>         potential_type;

  public:

    BondAngleInteraction() = default;
    ~BondAngleInteraction() = default;

    void
    calc_force(particle_type& p1, particle_type& p2, particle_type& p3,
               const potential_type& pot) const override;

    real_type
    calc_energy(const particle_type& p1, const particle_type& p2,
                const particle_type& p3, const potential_type& pot) const override;
};

// \frac{\partial\theta}{\partial\r_pre}
// = \frac{1}{|r_pre - r_mid| \sin\theta}
//   (- \frac{r_post - r_mid}{|r_post - r_mid|}
//    + \frac{r_pre  - r_mid}{|r_pre  - r_mid|}\cos\theta)
template<typename traitsT, typename boundaryT>
void BondAngleInteraction<traitsT, boundaryT>::calc_force(
        particle_type& p1, particle_type& p2, particle_type& p3,
        const potential_type& pot) const
{
    const coordinate_type r_ij =
        boundary_type::adjust_direction(p1.position - p2.position);
    const real_type       inv_len_r_ij = fast_inv_sqrt(length_sq(r_ij));
    const coordinate_type r_ij_reg     = r_ij * inv_len_r_ij;

    const coordinate_type r_kj =
        boundary_type::adjust_direction(p3.position - p2.position);
    const real_type       inv_len_r_kj = fast_inv_sqrt(length_sq(r_kj));
    const coordinate_type r_kj_reg     = r_kj * inv_len_r_kj;

    const real_type dot_ijk = dot_product(r_ij_reg, r_kj_reg);
    const real_type cos_theta =
        (-1. <= dot_ijk && dot_ijk <= 1.) ? dot_ijk : std::copysign(1.0, dot_ijk);

    const real_type theta = std::acos(cos_theta);

    const real_type coef = -pot.derivative(theta);

    const real_type sin_theta = std::sin(theta);
    const real_type coef_inv_sin = (sin_theta > constants<traits_type>::tolerance)
        ? coef / sin_theta : coef / constants<traits_type>::tolerance;

    const coordinate_type Fi =
        (coef_inv_sin * inv_len_r_ij) * (cos_theta * r_ij_reg - r_kj_reg);

    const coordinate_type Fk =
        (coef_inv_sin * inv_len_r_kj) * (cos_theta * r_kj_reg - r_ij_reg);

    p1.force += Fi;
    p2.force -= (Fi + Fk);
    p3.force += Fk;
    return;
}

template<typename traitsT, typename boundaryT>
typename BondAngleInteraction<traitsT, boundaryT>::real_type
BondAngleInteraction<traitsT, boundaryT>::calc_energy(
    const particle_type& p1, const particle_type& p2, const particle_type& p3,
    const potential_type& pot) const
{
    const coordinate_type v_2to1 =
        boundary_type::adjust_direction(p1.position - p2.position);
    const coordinate_type v_2to3 =
        boundary_type::adjust_direction(p3.position - p2.position);

    const real_type lensq_v21 = length_sq(v_2to1);
    const real_type lensq_v23 = length_sq(v_2to3);
    const real_type dot_v21_v23 = dot_product(v_2to1, v_2to3);

    const real_type dot_ijk = dot_v21_v23 * fast_inv_sqrt(lensq_v21 * lensq_v23);
    const real_type cos_theta = (-1. <= dot_ijk && dot_ijk <= 1.)
                                ? dot_ijk : std::copysign(1.0, dot_ijk);
    const real_type theta = std::acos(cos_theta);

    return pot.potential(theta);
}

}// mjolnir

#endif /* MJOLNIR_BOND_ANGL_INTERACTION */
