#ifndef MJOLNIR_CORE_BOND_ANGLE_INTERACTION
#define MJOLNIR_CORE_BOND_ANGLE_INTERACTION
#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/logger.hpp>
#include <cmath>

namespace mjolnir
{

template<typename traitsT, typename potentialT>
class BondAngleInteraction : public LocalInteractionBase<traitsT>
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

    typedef std::array<std::size_t, 3>          indices_type;
    typedef std::pair<indices_type, potentialT> potential_index_pair;
    typedef std::vector<potential_index_pair>   container_type;
    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;

    typedef math::constants<real_type> constant;

  public:

    BondAngleInteraction(const connection_kind_type kind,
                         const container_type& pot)
        : kind_(kind), potentials(pot)
    {}
    BondAngleInteraction(const connection_kind_type kind,
                         container_type&& pot)
        : kind_(kind), potentials(std::move(pot))
    {}
    ~BondAngleInteraction() override = default;

    void      calc_force (system_type&)        const noexcept override;
    real_type calc_energy(const system_type& ) const noexcept override;

    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_SCOPE(BondAngleInteraction::initialize(), 0);
        MJOLNIR_LOG_INFO("potential = ", potential_type::name(),
                         ", number of angles = ", potentials.size());
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
    {return "BondAngle:"_s + potential_type::name();}

    void write_topology(topology_type&) const override;

  private:
    connection_kind_type kind_;
    container_type potentials;
};

template<typename traitsT, typename potentialT>
void
BondAngleInteraction<traitsT, potentialT>::calc_force(system_type& sys) const noexcept
{
    for(const auto& idxp : this->potentials)
    {
        const std::size_t idx0 = idxp.first[0];
        const std::size_t idx1 = idxp.first[1];
        const std::size_t idx2 = idxp.first[2];

        const coordinate_type r_ij =
            sys.adjust_direction(sys[idx0].position - sys[idx1].position);

        const real_type       inv_len_r_ij = rsqrt(length_sq(r_ij));
        const coordinate_type r_ij_reg     = r_ij * inv_len_r_ij;

        const coordinate_type r_kj =
            sys.adjust_direction(sys[idx2].position - sys[idx1].position);

        const real_type       inv_len_r_kj = rsqrt(length_sq(r_kj));
        const coordinate_type r_kj_reg     = r_kj * inv_len_r_kj;

        const real_type dot_ijk   = dot_product(r_ij_reg, r_kj_reg);
        const real_type cos_theta = clamp(dot_ijk, real_type(-1.0), real_type(1.0));

        const real_type theta = std::acos(cos_theta);
        const real_type coef  = -(idxp.second.derivative(theta));

        const real_type sin_theta    = std::sin(theta);
        const real_type coef_inv_sin = (sin_theta > constant::tolerance) ?
                             coef / sin_theta : coef / constant::tolerance;

        const coordinate_type Fi =
            (coef_inv_sin * inv_len_r_ij) * (cos_theta * r_ij_reg - r_kj_reg);

        const coordinate_type Fk =
            (coef_inv_sin * inv_len_r_kj) * (cos_theta * r_kj_reg - r_ij_reg);

        sys[idx0].force += Fi;
        sys[idx1].force -= (Fi + Fk);
        sys[idx2].force += Fk;
    }
    return;
}

template<typename traitsT, typename potentialT>
typename BondAngleInteraction<traitsT, potentialT>::real_type
BondAngleInteraction<traitsT, potentialT>::calc_energy(
        const system_type& sys) const noexcept
{
    real_type E = 0.0;
    for(const auto& idxp : this->potentials)
    {
        const std::size_t idx0 = idxp.first[0];
        const std::size_t idx1 = idxp.first[1];
        const std::size_t idx2 = idxp.first[2];

        const coordinate_type v_2to1 =
            sys.adjust_direction(sys[idx0].position - sys[idx1].position);
        const coordinate_type v_2to3 =
            sys.adjust_direction(sys[idx2].position - sys[idx1].position);

        const real_type lensq_v21   = length_sq(v_2to1);
        const real_type lensq_v23   = length_sq(v_2to3);
        const real_type dot_v21_v23 = dot_product(v_2to1, v_2to3);

        const real_type dot_ijk   = dot_v21_v23 * rsqrt(lensq_v21 * lensq_v23);
        const real_type cos_theta = clamp(dot_ijk, real_type(-1.0), real_type(1.0));
        const real_type theta     = std::acos(cos_theta);

        E += idxp.second.potential(theta);
    }
    return E;
}

template<typename traitsT, typename potentialT>
void BondAngleInteraction<traitsT, potentialT>::write_topology(
        topology_type& topol) const
{
    if(this->kind_.empty() || this->kind_ == "none") {return;}

    for(const auto& idxp : this->potentials)
    {
        const auto i = idxp.first[0];
        const auto j = idxp.first[1];
        const auto k = idxp.first[2];
        topol.add_connection(i, j, this->kind_);
        topol.add_connection(j, k, this->kind_);
    }
    return;
}

}// mjolnir
#endif /* MJOLNIR_BOND_ANGL_INTERACTION */
