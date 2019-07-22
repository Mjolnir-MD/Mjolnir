#ifndef MJOLNIR_INTERACTION_BOND_ANGLE_INTERACTION_HPP
#define MJOLNIR_INTERACTION_BOND_ANGLE_INTERACTION_HPP
#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/logger.hpp>
#include <cmath>

namespace mjolnir
{

// calculate energy and force of bond angle interaction.
// The implementation is based on the book written by A. Ueda (2003) Shokabo.
template<typename traitsT, typename potentialT>
class BondAngleInteraction final : public LocalInteractionBase<traitsT>
{
  public:
    using traits_type          = traitsT;
    using potential_type       = potentialT;
    using base_type            = LocalInteractionBase<traits_type>;
    using real_type            = typename base_type::real_type;
    using coordinate_type      = typename base_type::coordinate_type;
    using system_type          = typename base_type::system_type;
    using topology_type        = typename base_type::topology_type;
    using connection_kind_type = typename base_type::connection_kind_type;

    using indices_type         = std::array<std::size_t, 3>;
    using potential_index_pair = std::pair<indices_type, potentialT>;
    using container_type       = std::vector<potential_index_pair>;
    using iterator             = typename container_type::iterator;
    using const_iterator       = typename container_type::const_iterator;

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

    void initialize(const system_type&) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
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

    // do nothing. this is used to reduce margin of neighbor list, and added
    // to this class for the consistency.
    void update_margin(const real_type, const system_type&) override {return;}

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
            sys.adjust_direction(sys.position(idx0) - sys.position(idx1));

        const real_type       inv_len_r_ij = math::rlength(r_ij);
        const coordinate_type r_ij_reg     = r_ij * inv_len_r_ij;

        const coordinate_type r_kj =
            sys.adjust_direction(sys.position(idx2) - sys.position(idx1));

        const real_type       inv_len_r_kj = math::rlength(r_kj);
        const coordinate_type r_kj_reg     = r_kj * inv_len_r_kj;

        const real_type dot_ijk   = math::dot_product(r_ij_reg, r_kj_reg);
        const real_type cos_theta = math::clamp(dot_ijk, real_type(-1.0), real_type(1.0));

        const real_type theta = std::acos(cos_theta);
        const real_type coef  = -(idxp.second.derivative(theta));

        const real_type sin_theta    = std::sin(theta);
        const real_type coef_inv_sin = (sin_theta > math::tolerance<real_type>()) ?
                             coef / sin_theta : coef / math::tolerance<real_type>();

        const coordinate_type Fi =
            (coef_inv_sin * inv_len_r_ij) * (cos_theta * r_ij_reg - r_kj_reg);

        const coordinate_type Fk =
            (coef_inv_sin * inv_len_r_kj) * (cos_theta * r_kj_reg - r_ij_reg);

        sys.force(idx0) += Fi;
        sys.force(idx1) -= (Fi + Fk);
        sys.force(idx2) += Fk;
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
            sys.adjust_direction(sys.position(idx0) - sys.position(idx1));
        const coordinate_type v_2to3 =
            sys.adjust_direction(sys.position(idx2) - sys.position(idx1));

        const real_type lensq_v21   = math::length_sq(v_2to1);
        const real_type lensq_v23   = math::length_sq(v_2to3);
        const real_type dot_v21_v23 = math::dot_product(v_2to1, v_2to3);

        const real_type dot_ijk   = dot_v21_v23 * math::rsqrt(lensq_v21 * lensq_v23);
        const real_type cos_theta = math::clamp(dot_ijk, real_type(-1.0), real_type(1.0));
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
        topol.add_connection(i, k, this->kind_);
        topol.add_connection(j, k, this->kind_);
    }
    return;
}

}// mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize BondAngleInteraction with LocalPotentials
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/potential/local/HarmonicPotential.hpp>
#include <mjolnir/potential/local/GaussianPotential.hpp>
#include <mjolnir/potential/local/FlexibleLocalAnglePotential.hpp>

namespace mjolnir
{

// harmonic
extern template class BondAngleInteraction<SimulatorTraits<double, UnlimitedBoundary>, HarmonicPotential<double>>;
extern template class BondAngleInteraction<SimulatorTraits<float,  UnlimitedBoundary>, HarmonicPotential<float> >;
extern template class BondAngleInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, HarmonicPotential<double>>;
extern template class BondAngleInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, HarmonicPotential<float> >;

// gaussian
extern template class BondAngleInteraction<SimulatorTraits<double, UnlimitedBoundary>, GaussianPotential<double>>;
extern template class BondAngleInteraction<SimulatorTraits<float,  UnlimitedBoundary>, GaussianPotential<float> >;
extern template class BondAngleInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>>;
extern template class BondAngleInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> >;

// FLP angle
extern template class BondAngleInteraction<SimulatorTraits<double, UnlimitedBoundary>, FlexibleLocalAnglePotential<double>>;
extern template class BondAngleInteraction<SimulatorTraits<float,  UnlimitedBoundary>, FlexibleLocalAnglePotential<float> >;
extern template class BondAngleInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, FlexibleLocalAnglePotential<double>>;
extern template class BondAngleInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, FlexibleLocalAnglePotential<float> >;

} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_BOND_ANGL_INTERACTION */
