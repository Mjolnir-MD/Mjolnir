#ifndef JARNGREIPR_GEOMETRY_DIHEDRAL
#define JARNGREIPR_GEOMETRY_DIHEDRAL
#include "angle.hpp"

namespace jarngreipr
{

template<typename coordT>
typename scalar_type_of<coordT>::type dihedral_angle(
        const coordT& p1, const coordT& p2, const coordT& p3, const coordT& p4)
{
    const auto r_ij = p1 - p2;
    const auto r_kj = p3 - p2;
    const auto r_lk = p4 - p3;
    const auto r_kj_lensq_inv = 1. / mjolnir::length_sq(r_kj);

    const auto n = mjolnir::cross_product(r_kj, -1e0 * r_lk);

    const auto R = r_ij - (mjolnir::dot_product(r_ij, r_kj) * r_kj_lensq_inv) * r_kj;
    const auto S = r_lk - (mjolnir::dot_product(r_lk, r_kj) * r_kj_lensq_inv) * r_kj;

    const auto cos_RS = mjolnir::dot_product(R, S) *
        mjolnir::fast_inv_sqrt(mjolnir::length_sq(R) * mjolnir::length_sq(S));

    const auto cos_phi = (-1e0 <= cos_RS && cos_RS <= 1e0)
                         ? cos_RS : std::copysign(1.0, cos_RS);
    const auto phi = std::copysign(std::acos(cos_phi), mjolnir::dot_product(r_ij, n));

    return cos_phi;
}




} // jarngreipr
#endif /* JARNGREIPR_GEOMETRY_DIHEDRAL */
