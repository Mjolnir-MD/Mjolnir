#ifndef MJOLNIR_FORCEFIELD_STOICHIOMETRIC_GLOBAL_STOICHIOMETRIC_INTERACTION_HPP
#define MJOLNIR_FORCEFIELD_STOICHIOMETRIC_GLOBAL_STOICHIOMETRIC_INTERACTION_HPP
#include <mjolnir/forcefield/stoichiometric/GlobalStoichiometricInteractionPotential.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/GlobalInteractionBase.hpp>
#include <mjolnir/core/SpatialPartitionBase.hpp>
#include <mjolnir/util/logger.hpp>
#include <memory>


namespace mjolnir
{

// The interaction which conserve stoichiometry.

template<typename traitsT>
class GlobalStoichiometricInteraction final : public GlobalInteractionBase<traitsT>
{
  public:

    using traits_type    = traitsT;
    using base_type      = GlobalInteractionBase<traits_type>;
    using real_type      = typename base_type::real_type;
    using system_type    = typename base_type::system_type;
    using topology_type  = typename base_type::topology_type;
    using potential_type = GlobalStoichiometricInteractionPotential<traits_type>;
    using partition_type = SpatialPartition<traitsT, potential_type>;

  public:
    GlobalStoichiometricInteraction(potential_type&& pot, partition_type&& part)
        : potential_(std::move(pot)), partition_(std::move(part))
    {}
    ~GlobalStoichiometricInteraction() {}

    void initialize(const system_type& sys, const topology_type& topol) override
    {
    }

    void update(const system_type&, const topology_type&) override
    {
    }

    void reduce_margin(const real_type, const system_type&) override
    {
    }
    void scale_margin(const real_type, const system_type&) override
    {
    }

    void      calc_force (system_type&)       const noexcept override;
    real_type calc_energy(const system_type&) const noexcept override;

    std::string name() const override {return "Stoichiometric";}

    potential_type const& potential() const noexcept {return potential_;}
    potential_type &      potential()       noexcept {return potential_;}

    base_type* clone() const override
    {
        return new GlobalStoichiometricInteraction(
                potential_type(potential_), partition_type(partition_));
    }

  private:

    potential_type potential_;
    partition_type partition_;
};

template<typename traitsT>
void GlobalStoichiometricInteraction<traitsT>::calc_force(
        system_type& sys) const noexcept
{
    return ;
}

template<typename traitsT>
typename GlobalStoichiometricInteraction<traitsT>::real_type
GlobalStoichiometricInteraction<traitsT>::calc_energy(
        const system_type& sys) const noexcept
{
    return 0.0;
}
} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize major use-cases
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

namespace mjolnir
{

extern template class GlobalStoichiometricInteraction<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class GlobalStoichiometricInteraction<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class GlobalStoichiometricInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class GlobalStoichiometricInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_FORCEFIELD_STOICHIOMETRIC_GLOBAL_STOICHIOMETRIC_INTERACTION_HPP */
