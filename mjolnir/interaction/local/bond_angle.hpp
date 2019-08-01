#ifndef MJOLNIR_INTERACTION_LOCAL_BOND_ANGLE_HPP
#define MJOLNIR_INTERACTION_LOCAL_BOND_ANGLE_HPP
#include <mjolnir/math/math.hpp>

namespace mjolnir
{

// it calculates bond angle and directions of forces.
//
//        o r0            |
//       /                |
//   r1 o) <- bond angle  |
//       \                |
//        o r2            |
//
// Note for developers: At first, this function was splitted from BondAngle-
// Interaction. Then it turned out calling this from BondAngleInteraction::
// calc_force slows down the simulation in ~5% by a micro benchmark using
// AICG2+ protein forcefield. It was a surprising result.
//
// Also I found that, unlike dihedral_angle case, `attribute((always_inline))`
// DOES NOT recover the efficiency. I have no idea why.
//
// Eventually I decided to keep the BondAngleInteraction and implement
// bond angle calculater here independently. The reason is that
// - It takes a lot of efforts to reduce the run time 5% by other tricks.
// - We can test this function also.
// - Splitting this from BondAngleInteraction introduces an additional
//   complexity to the BondAngleInteraction because other LocalInteraction
//   classes (currently) does not have the corresponding function.
template<typename traitsT>
std::pair<typename traitsT::real_type,
          std::array<typename traitsT::coordinate_type, 3>>
calc_bond_angle_force(const System<traitsT>& sys,
        const typename traitsT::coordinate_type& r0,
        const typename traitsT::coordinate_type& r1,
        const typename traitsT::coordinate_type& r2) noexcept
{
    using real_type = typename traitsT::real_type;
    using coordinate_type = typename traitsT::coordinate_type;
    const auto r_ij = sys.adjust_direction(r0 - r1);
    const auto r_kj = sys.adjust_direction(r2 - r1);

    const auto inv_len_r_ij = math::rlength(r_ij);
    const auto r_ij_reg     = r_ij * inv_len_r_ij;

    const auto inv_len_r_kj = math::rlength(r_kj);
    const auto r_kj_reg     = r_kj * inv_len_r_kj;

    const auto dot_ijk   = math::dot_product(r_ij_reg, r_kj_reg);
    const auto cos_theta = math::clamp<real_type>(dot_ijk, -1, 1);

    const auto theta = std::acos(cos_theta);

    constexpr auto tol   = math::abs_tolerance<real_type>();
    const auto sin_theta = std::sin(theta);
    const auto inv_sin   = (sin_theta > tol) ? real_type(1.0) / sin_theta :
                                               real_type(1.0) / tol;
    std::array<coordinate_type, 3> F;
    // take care about the indices.
    F[0] = (inv_sin * inv_len_r_ij) * (cos_theta * r_ij_reg - r_kj_reg);
    F[2] = (inv_sin * inv_len_r_kj) * (cos_theta * r_kj_reg - r_ij_reg);
    F[1] = -1 * (F[0] + F[2]);

    return std::make_pair(theta, F);
}

template<typename traitsT>
typename traitsT::real_type
calc_bond_angle(const System<traitsT>& sys,
        const typename traitsT::coordinate_type& r0,
        const typename traitsT::coordinate_type& r1,
        const typename traitsT::coordinate_type& r2) noexcept
{
    using real_type = typename traitsT::real_type;

    const auto r_ij = sys.adjust_direction(r0 - r1);
    const auto r_kj = sys.adjust_direction(r2 - r1);

    const auto dot_ijk   = math::dot_product(r_ij, r_kj) *
                           math::rlength(r_ij) * math::rlength(r_kj);
    const auto cos_theta = math::clamp<real_type>(dot_ijk, -1, 1);

    return std::acos(cos_theta);
}

} // mjolnir
#endif// MJOLNIR_INTERACTION_LOCAL_BOND_ANGLE_HPP
