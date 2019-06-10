#ifndef MJOLNIR_INTERACTION_LOCAL_DIHEDRAL_ANGLE_HPP
#define MJOLNIR_INTERACTION_LOCAL_DIHEDRAL_ANGLE_HPP
#include <mjolnir/math/math.hpp>

namespace mjolnir
{

// it calculates dihedral angle and direction of forces.
//
//        o r0
//       /
//   r1 o
//      (\) <- dihedral angle
//        o r2
//       /
//   r3 o
//
// Note for developers: At first, this function was splitted from Dihedral-
// AngleInteraction. Then it turned out calling this from DihedralAngle-
// Interaction::calc_force slows down the simulation in ~5% by a micro
// benchmark using AICG2+ protein forcefield. It was a surprising result.
// Also I found that `__attribute__((always_inline))` recovers the efficiency,
// while `inline` does not. Possible reasons include
// - The function was inlined before other optimization procedures and it makes
//   the successive optimization effective.
// - Simply `inline` was ignored by a compiler because of the size.
// But if the reason was the optimization flow, it means that this recovery
// depends on the version of the compiler.
//
// Eventually I decided to keep the DihedralAngleInteraction and implement
// dihedral angle calculater here independently. The reason is that
// - It takes a lot of efforts to reduce the run time 5% by other tricks.
// - We can test this function also.
// - Splitting this from DihedralAngleInteraction introduces an additional
//   complexity to the DihedralAngleInteraction because other LocalInteraction
//   classes (currently) does not have the corresponding function.
template<typename traitsT>
std::pair<typename traitsT::real_type,
          std::array<typename traitsT::coordinate_type, 4>>
calc_dihedral_angle_force(const System<traitsT>& sys,
        const typename traitsT::coordinate_type& r0,
        const typename traitsT::coordinate_type& r1,
        const typename traitsT::coordinate_type& r2,
        const typename traitsT::coordinate_type& r3) noexcept
{
    using real_type = typename traitsT::real_type;
    const auto r_ij = sys.adjust_direction(r0 - r1);
    const auto r_kj = sys.adjust_direction(r2 - r1);
    const auto r_lk = sys.adjust_direction(r3 - r2);
    const auto r_kl = -1 * r_lk;

    const auto r_kj_lensq  = math::length_sq(r_kj);
    const auto r_kj_rlen   = math::rsqrt(r_kj_lensq);
    const auto r_kj_rlensq = r_kj_rlen * r_kj_rlen;
    const auto r_kj_len    = r_kj_rlen * r_kj_lensq;

    const auto m       = math::cross_product(r_ij, r_kj);
    const auto n       = math::cross_product(r_kj, r_kl);
    const auto m_lensq = math::length_sq(m);
    const auto n_lensq = math::length_sq(n);

    const auto R = r_ij - (math::dot_product(r_ij, r_kj) * r_kj_rlensq) * r_kj;
    const auto S = r_lk - (math::dot_product(r_lk, r_kj) * r_kj_rlensq) * r_kj;
    const auto R_lensq = math::length_sq(R);
    const auto S_lensq = math::length_sq(S);

    const auto dot_RS  = math::dot_product(R, S) * math::rsqrt(R_lensq * S_lensq);
    const auto cos_phi = math::clamp<real_type>(dot_RS, -1, 1);
    const auto phi     = std::copysign(std::acos(cos_phi),
                                       math::dot_product(r_ij, n));
    using coordinate_type = typename traitsT::coordinate_type;
    std::array<coordinate_type, 4> f;

    f[0] = ( r_kj_len / m_lensq) * m;
    f[3] = (-r_kj_len / n_lensq) * n;

    const auto coef_ijk = math::dot_product(r_ij, r_kj) * r_kj_rlensq;
    const auto coef_jkl = math::dot_product(r_kl, r_kj) * r_kj_rlensq;

    f[1] = (coef_ijk - 1) * f[0] - coef_jkl * f[3];
    f[2] = (coef_jkl - 1) * f[3] - coef_ijk * f[0];

    return std::make_pair(phi, f);
}

template<typename traitsT>
typename traitsT::real_type
calc_dihedral_angle(const System<traitsT>& sys,
        const typename traitsT::coordinate_type& r0,
        const typename traitsT::coordinate_type& r1,
        const typename traitsT::coordinate_type& r2,
        const typename traitsT::coordinate_type& r3) noexcept
{
    using real_type = typename traitsT::real_type;
    const auto r_ij = sys.adjust_direction(r0 - r1);
    const auto r_kj = sys.adjust_direction(r2 - r1);
    const auto r_lk = sys.adjust_direction(r3 - r2);

    const auto r_kj_lensq  = math::length_sq(r_kj);
    const auto r_kj_rlen   = math::rsqrt(r_kj_lensq);
    const auto r_kj_rlensq = r_kj_rlen * r_kj_rlen;

    const auto n = math::cross_product(r_kj, -1 * r_lk);

    const auto R = r_ij - (math::dot_product(r_ij, r_kj) * r_kj_rlensq) * r_kj;
    const auto S = r_lk - (math::dot_product(r_lk, r_kj) * r_kj_rlensq) * r_kj;
    const auto R_lensq = math::length_sq(R);
    const auto S_lensq = math::length_sq(S);

    const auto dot_RS  = math::dot_product(R, S) * math::rsqrt(R_lensq * S_lensq);
    const auto cos_phi = math::clamp<real_type>(dot_RS, -1, 1);

    return std::copysign(std::acos(cos_phi), math::dot_product(r_ij, n));
}
} // mjolnir
#endif // MJOLNIR_INTERACTION_LOCAL_DIHEDRAL_ANGLE_HPP
