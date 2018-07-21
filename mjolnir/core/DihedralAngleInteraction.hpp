#ifndef MJOLNIR_CORE_DIHEDRAL_ANGLE_INTERACTION
#define MJOLNIR_CORE_DIHEDRAL_ANGLE_INTERACTION
#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/math/rsqrt.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/logger.hpp>
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
    typedef typename base_type::topology_type   topology_type;
    typedef typename base_type::connection_kind_type connection_kind_type;

    typedef std::array<std::size_t, 4>          indices_type;
    typedef std::pair<indices_type, potentialT> potential_index_pair;
    typedef std::vector<potential_index_pair>   container_type;
    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;

  public:

    DihedralAngleInteraction(const connection_kind_type kind,
                             const container_type& pot)
        : kind_(kind), potentials(pot)
    {}
    DihedralAngleInteraction(const connection_kind_type kind,
                             container_type&& pot)
        : kind_(kind), potentials(std::move(pot))
    {}
    ~DihedralAngleInteraction() = default;

    void      calc_force (system_type&)       const noexcept override;
    real_type calc_energy(const system_type&) const noexcept override;

    void initialize(const system_type& sys, const real_type dt) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_SCOPE(DihedralAngleInteraction::initialize(), 0);
        MJOLNIR_LOG_INFO("potential = ", potential_type::name(),
                         ", N = ", potentials.size());
        return;
    }

    void update(const system_type& sys) override
    {
        for(auto& item : potentials)
        {
            item.second.update(sys);
        }
    }

    std::string name() const override
    {return "DihedralAngle:"_str + potential_type::name();}

    void write_topology(topology_type&) const override;

   private:
    connection_kind_type kind_;
    container_type potentials;
};


template<typename traitsT, typename pT>
void
DihedralAngleInteraction<traitsT, pT>::calc_force(system_type& sys) const noexcept
{
    for(const auto& idxp : this->potentials)
    {
        const std::size_t idx0 = idxp.first[0];
        const std::size_t idx1 = idxp.first[1];
        const std::size_t idx2 = idxp.first[2];
        const std::size_t idx3 = idxp.first[3];

        const coordinate_type r_ij =
            sys.adjust_direction(sys[idx0].position - sys[idx1].position);
        const coordinate_type r_kj =
            sys.adjust_direction(sys[idx2].position - sys[idx1].position);
        const coordinate_type r_lk =
            sys.adjust_direction(sys[idx3].position - sys[idx2].position);
        const coordinate_type r_kl = -1e0 * r_lk;

        const real_type r_kj_lensq  = length_sq(r_kj);
        const real_type r_kj_rlen   = rsqrt(r_kj_lensq);
        const real_type r_kj_rlensq = r_kj_rlen * r_kj_rlen;
        const real_type r_kj_len    = r_kj_rlen * r_kj_lensq;

        const coordinate_type m = cross_product(r_ij, r_kj);
        const coordinate_type n = cross_product(r_kj, r_kl);
        const real_type m_lensq = length_sq(m);
        const real_type n_lensq = length_sq(n);

        const coordinate_type R =
            r_ij - (dot_product(r_ij, r_kj) * r_kj_rlensq) * r_kj;
        const coordinate_type S =
            r_lk - (dot_product(r_lk, r_kj) * r_kj_rlensq) * r_kj;

        const real_type R_lensq = length_sq(R);
        const real_type S_lensq = length_sq(S);

        const real_type dot_RS  = dot_product(R, S) * rsqrt(R_lensq * S_lensq);
        const real_type cos_phi = (-1. <= dot_RS && dot_RS <= 1.)
                                  ? dot_RS : std::copysign(1.0, dot_RS);
        const real_type phi =
            std::copysign(std::acos(cos_phi), dot_product(r_ij, n));

        // -dV / dphi
        const real_type coef = -(idxp.second.derivative(phi));

        const coordinate_type Fi = ( coef * r_kj_len / m_lensq) * m;
        const coordinate_type Fl = (-coef * r_kj_len / n_lensq) * n;

        const real_type coef_ijk = dot_product(r_ij, r_kj) * r_kj_rlensq;
        const real_type coef_jkl = dot_product(r_kl, r_kj) * r_kj_rlensq;

        sys[idx0].force += Fi;
        sys[idx1].force += (coef_ijk - 1.0) * Fi - coef_jkl * Fl;
        sys[idx2].force += (coef_jkl - 1.0) * Fl - coef_ijk * Fi;
        sys[idx3].force += Fl;
    }
    return;
}

template<typename traitsT, typename potentialT>
typename DihedralAngleInteraction<traitsT, potentialT>::real_type
DihedralAngleInteraction<traitsT, potentialT>::calc_energy(
        const system_type& sys) const noexcept
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

        const real_type dot_RS = dot_product(R, S) * rsqrt(R_lensq * S_lensq);
        const real_type cos_phi = (-1e0 <= dot_RS && dot_RS <= 1e0)
                                  ? dot_RS : std::copysign(1.0, dot_RS);
        const real_type phi =
            std::copysign(std::acos(cos_phi), dot_product(r_ij, n));

        E += idxp.second.potential(phi);
    }
    return E;
}

template<typename traitsT, typename potentialT>
void DihedralAngleInteraction<traitsT, potentialT>::write_topology(
        topology_type& topol) const
{
    if(this->kind_.empty() || this->kind_ == "none") {return;}

    for(const auto& idxp : this->potentials)
    {
        const auto i = idxp.first[0];
        const auto j = idxp.first[1];
        const auto k = idxp.first[2];
        const auto l = idxp.first[3];
        topol.add_connection(i, j, this->kind_);
        topol.add_connection(j, k, this->kind_);
        topol.add_connection(k, l, this->kind_);
    }
    return;
}

}// mjolnir
#endif /* MJOLNIR_DIHEDRAL_ANGLE_INTERACTION */
