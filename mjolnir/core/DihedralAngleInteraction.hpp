#ifndef MJOLNIR_DIHEDRAL_ANGLE_INTERACTION
#define MJOLNIR_DIHEDRAL_ANGLE_INTERACTION
#include "LocalInteractionBase.hpp"
#include "BoundaryCondition.hpp"
#include <mjolnir/math/fast_inv_sqrt.hpp>
#include <cmath>

namespace mjolnir
{

/*! @brief calculate energy and force of Bond length type local interaction */
template<typename traitsT, typename potentialT,
         typename boundaryT = UnlimitedBoundary<traitsT>>
class DihedralAngleInteraction : public LocalInteractionBase<traitsT, 4>
{
  public:
    typedef traitsT    traits_type;
    typedef potentialT potential_type;
    typedef boundaryT  boundary_type;
    typedef LocalInteractionBase<traits_type, 4> base_type;
    typedef typename base_type::time_type        time_type;
    typedef typename base_type::real_type        real_type;
    typedef typename base_type::coordinate_type  coordinate_type;
    typedef typename base_type::particle_type    particle_type;

  public:

    DihedralAngleInteraction() = default;
    DihedralAngleInteraction(const potential_type& pot): potential(pot){}
    DihedralAngleInteraction(potential_type&& pot)
        : potential(std::forward<potential_type>(pot)){}
    ~DihedralAngleInteraction() = default;
    DihedralAngleInteraction(const DihedralAngleInteraction&) = default;
    DihedralAngleInteraction(DihedralAngleInteraction&&) = default;
    DihedralAngleInteraction& operator=(const DihedralAngleInteraction&) = default;
    DihedralAngleInteraction& operator=(DihedralAngleInteraction&&) = default;

    void
    calc_force(particle_type& p1, particle_type& p2, particle_type& p3,
               particle_type& p4) const override;

    real_type
    calc_energy(const particle_type& p1, const particle_type& p2,
                const particle_type& p3, const particle_type& p4) const override;
  private:
    potential_type potential;
};

template<typename traitsT, typename potentialT, typename boundaryT>
void DihedralAngleInteraction<traitsT, potentialT, boundaryT>::calc_force(
    particle_type& p1, particle_type& p2, particle_type& p3, particle_type& p4) const
{
    const coordinate_type r_ij =
        boundary_type::adjust_direction(p1.position - p2.position);
    const coordinate_type r_kj =
        boundary_type::adjust_direction(p3.position - p2.position);
    const coordinate_type r_lk =
        boundary_type::adjust_direction(p4.position - p3.position);
    const coordinate_type r_kl = -1e0 * r_lk;

    const real_type r_kj_lensq     = length_sq(r_kj);
    const real_type r_kj_lensq_inv = 1. / r_kj_lensq;
    const real_type r_kj_len       = std::sqrt(r_kj_lensq);

    const coordinate_type m = cross_product(r_ij, r_kj);
    const coordinate_type n = cross_product(r_kj, r_kl);
    const real_type m_lensq = length_sq(m);
    const real_type n_lensq = length_sq(n);

    const coordinate_type R = r_ij -
                              (dot_product(r_ij, r_kj) * r_kj_lensq_inv) * r_kj;
    const coordinate_type S = r_lk -
                              (dot_product(r_lk, r_kj) * r_kj_lensq_inv) * r_kj;
    const real_type R_lensq = length_sq(R);
    const real_type S_lensq = length_sq(S);

    const real_type dot_RS  = dot_product(R, S) *
                              fast_inv_sqrt(R_lensq * S_lensq);
    const real_type cos_phi = (-1. <= dot_RS && dot_RS <= 1.)
                              ? dot_RS : std::copysign(1.0, dot_RS);
    const real_type phi = std::copysign(std::acos(cos_phi), dot_product(r_ij, n));

    // -dV / dphi
    const real_type coef = -potential.derivative(phi);

    const coordinate_type Fi = ( coef * r_kj_len / m_lensq) * m;
    const coordinate_type Fl = (-coef * r_kj_len / n_lensq) * n;

    const real_type coef_ijk = dot_product(r_ij, r_kj) * r_kj_lensq_inv;
    const real_type coef_jkl = dot_product(r_kl, r_kj) * r_kj_lensq_inv;

    p1.force += Fi;
    p2.force += (coef_ijk - 1e0) * Fi - coef_jkl * Fl;
    p3.force += (coef_jkl - 1e0) * Fl - coef_ijk * Fi;
    p4.force += Fl;
    return;
}

template<typename traitsT, typename potentialT, typename boundaryT>
typename DihedralAngleInteraction<traitsT, potentialT, boundaryT>::real_type
DihedralAngleInteraction<traitsT, potentialT, boundaryT>::calc_energy(
        const particle_type& p1, const particle_type& p2,
        const particle_type& p3, const particle_type& p4) const
{
    const coordinate_type r_ij =
        boundary_type::adjust_direction(p1.position - p2.position);
    const coordinate_type r_kj =
        boundary_type::adjust_direction(p3.position - p2.position);
    const coordinate_type r_lk =
        boundary_type::adjust_direction(p4.position - p3.position);
    const real_type r_kj_lensq_inv = 1. / length_sq(r_kj);

    const coordinate_type n = cross_product(r_kj, -1e0 * r_lk);

    const coordinate_type R = r_ij -
                              (dot_product(r_ij, r_kj) * r_kj_lensq_inv) * r_kj;
    const coordinate_type S = r_lk -
                              (dot_product(r_lk, r_kj) * r_kj_lensq_inv) * r_kj;
    const real_type R_lensq = length_sq(R);
    const real_type S_lensq = length_sq(S);

    const real_type dot_RS = dot_product(R, S) *
                             fast_inv_sqrt(R_lensq * S_lensq);
    const real_type cos_phi = (-1e0 <= dot_RS && dot_RS <= 1e0)
                              ? dot_RS : std::copysign(1.0, dot_RS);
    const real_type phi = std::copysign(std::acos(cos_phi), dot_product(r_ij, n));

    return potential.potential(phi);
}

}// mjolnir
#endif /* MJOLNIR_DIHEDRAL_ANGLE_INTERACTION */
