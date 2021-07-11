#ifndef MJOLNIR_INTERACTION_BOND_LENGTH_INTERACTION_HPP
#define MJOLNIR_INTERACTION_BOND_LENGTH_INTERACTION_HPP
#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/logger.hpp>
#include <cmath>
#include <iostream>

namespace mjolnir
{

// Bond Length dependent interaction.
//
// For the following potentials, there are specialized definitions to speedup
// calculation.
// - GoContactPotential
//   - BondLengthGoContactInteraction.hpp
// By combining potential calculations, we can omit `sqrt()` that is usually
// used to calculate the distance between particles.
//
template<typename traitsT, typename potentialT>
class BondLengthInteraction final : public LocalInteractionBase<traitsT>
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

    using indices_type         = std::array<std::size_t, 2>;
    using potential_index_pair = std::pair<indices_type, potentialT>;
    using container_type       = std::vector<potential_index_pair>;
    using iterator             = typename container_type::iterator;
    using const_iterator       = typename container_type::const_iterator;

  public:

    BondLengthInteraction(const connection_kind_type kind,
                          const container_type& pot)
        : kind_(kind), potentials_(pot)
    {}
    BondLengthInteraction(const connection_kind_type kind,
                          container_type&& pot)
        : kind_(kind), potentials_(std::move(pot))
    {}
    ~BondLengthInteraction() override {}

    real_type calc_energy(const system_type&)     const noexcept override;

    void      calc_force           (system_type&) const noexcept override;
    void      calc_force_and_virial(system_type&) const noexcept override;
    real_type calc_force_and_energy(system_type&) const noexcept override;

    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential = ", potential_type::name(),
                         ", number of bonds = ", potentials_.size());
        for(auto& potential : this->potentials_)
        {
            potential.second.initialize(sys);
        }
        return;
    }

    void update(const system_type& sys) override
    {
        for(auto& item : potentials_)
        {
            item.second.update(sys);
        }
    }

    // do nothing. this is used to reduce margin of neighbor list, and added
    // to this class for the consistency.
    void reduce_margin(const real_type, const system_type&) override {return;}
    void  scale_margin(const real_type, const system_type&) override {return;}

    std::string name() const override
    {return "BondLength:"_s + potential_type::name();}

    void write_topology(topology_type&) const override;

    container_type const& potentials() const noexcept {return potentials_;}
    container_type&       potentials()       noexcept {return potentials_;}

    base_type* clone() const override
    {
        return new BondLengthInteraction(kind_, container_type(potentials_));
    }

  private:
    connection_kind_type kind_;
    container_type potentials_;

#ifdef MJOLNIR_WITH_OPENMP
    // OpenMP implementation uses its own implementation to run it in parallel.
    // So this implementation should not be instanciated with OpenMP Traits.
    static_assert(!is_openmp_simulator_traits<traits_type>::value,
                  "this is the default implementation, not for OpenMP");
#endif
};

template<typename traitsT, typename potentialT>
typename BondLengthInteraction<traitsT, potentialT>::real_type
BondLengthInteraction<traitsT, potentialT>::calc_energy(
        const system_type& sys) const noexcept
{
    real_type E = 0.;
    for(const auto& idxp : this->potentials_)
    {
        E += idxp.second.potential(math::length(sys.adjust_direction(
                sys.position(idxp.first[0]), sys.position(idxp.first[1]))));
    }
    return E;
}

template<typename traitsT, typename potentialT>
void BondLengthInteraction<traitsT, potentialT>::calc_force(
        system_type& sys) const noexcept
{
    for(const auto& idxp : this->potentials_)
    {
        const std::size_t idx0 = idxp.first[0];
        const std::size_t idx1 = idxp.first[1];

        const auto dpos = // from r0 -> r1 = r1 - r0
            sys.adjust_direction(sys.position(idx0), sys.position(idx1));

        const real_type len2 = math::length_sq(dpos); // l^2
        const real_type rlen = math::rsqrt(len2);     // 1/l
        const real_type force = -1 * idxp.second.derivative(len2 * rlen);
        // here, L^2 * (1 / L) = L.

        const coordinate_type f = dpos * (force * rlen);
        sys.force(idx0) -= f;
        sys.force(idx1) += f;
    }
    return;
}

template<typename traitsT, typename potentialT>
void BondLengthInteraction<traitsT, potentialT>::calc_force_and_virial(
        system_type& sys) const noexcept
{
    for(const auto& idxp : this->potentials_)
    {
        const std::size_t idx0 = idxp.first[0];
        const std::size_t idx1 = idxp.first[1];

        const auto dpos = // from r0 -> r1 = r1 - r0
            sys.adjust_direction(sys.position(idx0), sys.position(idx1));

        const real_type len2 = math::length_sq(dpos); // l^2
        const real_type rlen = math::rsqrt(len2);     // 1/l
        const real_type force = -1 * idxp.second.derivative(len2 * rlen);
        // here, L^2 * (1 / L) = L.

        const coordinate_type f = dpos * (force * rlen);
        sys.force(idx0) -= f;
        sys.force(idx1) += f;

        sys.virial() += math::tensor_product(dpos, f);
    }
    return;
}

template<typename traitsT, typename potentialT>
typename BondLengthInteraction<traitsT, potentialT>::real_type
BondLengthInteraction<traitsT, potentialT>::calc_force_and_energy(
        system_type& sys) const noexcept
{
    real_type energy = 0;
    for(const auto& idxp : this->potentials_)
    {
        const std::size_t idx0 = idxp.first[0];
        const std::size_t idx1 = idxp.first[1];

        const auto dpos =
            sys.adjust_direction(sys.position(idx0), sys.position(idx1));

        const real_type len2 = math::length_sq(dpos); // l^2
        const real_type rlen = math::rsqrt(len2);     // 1/l
        const real_type len  = len2 * rlen;
        const real_type force = -1 * idxp.second.derivative(len);
        energy += idxp.second.potential(len);

        const coordinate_type f = dpos * (force * rlen);
        sys.force(idx0) -= f;
        sys.force(idx1) += f;

        sys.virial() += math::tensor_product(dpos, f);
    }
    return energy;
}

template<typename traitsT, typename potentialT>
void BondLengthInteraction<traitsT, potentialT>::write_topology(
        topology_type& topol) const
{
    if(this->kind_.empty() || this->kind_ == "none") {return;}

    for(const auto& idxp : this->potentials_)
    {
        const auto i = idxp.first[0];
        const auto j = idxp.first[1];
        topol.add_connection(i, j, this->kind_);
    }
    return;
}

} // mjolnir

#include <mjolnir/forcefield/local/BondLengthGoContactInteraction.hpp>

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/local/HarmonicPotential.hpp>
#include <mjolnir/forcefield/local/GoContactPotential.hpp>
#include <mjolnir/forcefield/local/GaussianPotential.hpp>

namespace mjolnir
{

// harmonic
extern template class BondLengthInteraction<SimulatorTraits<double, UnlimitedBoundary       >, HarmonicPotential<double>>;
extern template class BondLengthInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, HarmonicPotential<float> >;
extern template class BondLengthInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, HarmonicPotential<double>>;
extern template class BondLengthInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, HarmonicPotential<float> >;

// The case of GoContactPotential is manually specialized in BondLengthGoContactInteraction.hpp.
// DO NOT specialize it here.

// gaussian
extern template class BondLengthInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>>;
extern template class BondLengthInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> >;
extern template class BondLengthInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>>;
extern template class BondLengthInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> >;

} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_BOND_LENGTH_INTERACTION */
