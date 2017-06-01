#ifndef MJOLNIR_DIHEDRAL_ANGLE_INTERACTION
#define MJOLNIR_DIHEDRAL_ANGLE_INTERACTION
#include "LocalInteractionBase.hpp"
#include <mjolnir/math/fast_inv_sqrt.hpp>
#include <cmath>

namespace mjolnir
{

/*! @brief calculate energy and force of Bond length type local interaction */
template<typename traitsT, typename potentialT>
class DihedralAngleInteraction : public LocalInteractionBase<traitsT>
{
  public:
    typedef traitsT    traits_type;
    typedef potentialT potential_type;
    typedef LocalInteractionBase<traits_type>   base_type;
    typedef typename base_type::real_type       real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::system_type     system_type;
    typedef typename base_type::particle_type   particle_type;
    typedef std::array<std::size_t, 4>          indices_type;
    typedef std::pair<indices_type, potentialT> potential_index_pair;
    typedef std::vector<potential_index_pair>   container_type;

  public:

    DihedralAngleInteraction() = default;
    ~DihedralAngleInteraction() = default;

    DihedralAngleInteraction(const container_type& pot): potentials(pot){}
    DihedralAngleInteraction(container_type&& pot): potentials(std::move(pot)){}

    void      calc_force (system_type&)       const override;
    real_type calc_energy(const system_type&) const override;

   private:
    container_type potentials;
};


template<typename traitsT, typename pT>
void
DihedralAngleInteraction<traitsT, pT>::calc_force(system_type& sys) const
{
    for(const auto& idxp : this->potentials)
    {
        const coordinate_type r_ij = sys.adjust_direction(
                sys[idxp.first[0]].position - sys[idxp.first[1]].position);
        const coordinate_type r_kj = sys.adjust_direction(
                sys[idxp.first[2]].position - sys[idxp.first[1]].position);
        const coordinate_type r_lk = sys.adjust_direction(
                sys[idxp.first[3]].position - sys[idxp.first[2]].position);
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
        const real_type phi = std::copysign(std::acos(cos_phi),
                                            dot_product(r_ij, n));

        // -dV / dphi
        const real_type coef = -(idxp.second.derivative(phi));

        const coordinate_type Fi = ( coef * r_kj_len / m_lensq) * m;
        const coordinate_type Fl = (-coef * r_kj_len / n_lensq) * n;

        const real_type coef_ijk = dot_product(r_ij, r_kj) * r_kj_lensq_inv;
        const real_type coef_jkl = dot_product(r_kl, r_kj) * r_kj_lensq_inv;

        sys[idxp.first[0]].force += Fi;
        sys[idxp.first[1]].force += (coef_ijk - 1e0) * Fi - coef_jkl * Fl;
        sys[idxp.first[2]].force += (coef_jkl - 1e0) * Fl - coef_ijk * Fi;
        sys[idxp.first[3]].force += Fl;
    }
    return;
}

template<typename traitsT, typename potentialT>
typename DihedralAngleInteraction<traitsT, potentialT>::real_type
DihedralAngleInteraction<traitsT, potentialT>::calc_energy(
        const system_type& sys) const
{
    real_type E = 0.;
    for(const auto& idxp : this->potentials)
    {
        const coordinate_type r_ij = sys.adjust_direction(
                sys[idxp.first[0]].position - sys[idxp.first[1]].position);
        const coordinate_type r_kj = sys.adjust_direction(
                sys[idxp.first[2]].position - sys[idxp.first[1]].position);
        const coordinate_type r_lk = sys.adjust_direction(
                sys[idxp.first[3]].position - sys[idxp.first[2]].position);
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
        const real_type phi =
            std::copysign(std::acos(cos_phi), dot_product(r_ij, n));

        E += idxp.second.potential(phi);
    }
    return E;
}

}// mjolnir
#endif /* MJOLNIR_DIHEDRAL_ANGLE_INTERACTION */
