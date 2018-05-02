#ifndef JARNGREIPR_GEOMETRY_DIHEDRAL
#define JARNGREIPR_GEOMETRY_DIHEDRAL
#include <jarngreipr/geometry/angle.hpp>

namespace jarngreipr
{

/*! @brief Calculate dihedral angle formed by p1, p2, p3, and p4.
 *      o p1   |            o p1
 *     /       |     dih _ /
 * p2 o        |        / /
 *     \       | p4 o----o p2 & p3
 *      o p3   |
 *     /       |
 * p4 o        |
 *             |
 * */
template<typename realT>
realT dihedral_angle(const mjolnir::Vector<realT, 3>& p1,
                     const mjolnir::Vector<realT, 3>& p2,
                     const mjolnir::Vector<realT, 3>& p3,
                     const mjolnir::Vector<realT, 3>& p4)
{
    const auto r_ij = p1 - p2;
    const auto r_kj = p3 - p2;
    const auto r_lk = p4 - p3;
    const auto r_kj_lensq_inv = 1. / mjolnir::length_sq(r_kj);

    const auto n = mjolnir::cross_product(r_kj, -1e0 * r_lk);

    const auto R = r_ij - (mjolnir::dot_product(r_ij, r_kj) * r_kj_lensq_inv) * r_kj;
    const auto S = r_lk - (mjolnir::dot_product(r_lk, r_kj) * r_kj_lensq_inv) * r_kj;

    const auto cos_RS = mjolnir::dot_product(R, S) * mjolnir::rsqrt(
            mjolnir::length_sq(R) * mjolnir::length_sq(S));

    const auto cos_phi = (-1e0 <= cos_RS && cos_RS <= 1e0)
                         ? cos_RS : std::copysign(1.0, cos_RS);

    return std::copysign(std::acos(cos_phi), dot_product(r_ij, n));
}

} // jarngreipr
#endif /* JARNGREIPR_GEOMETRY_DIHEDRAL */
