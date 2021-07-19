#ifndef MJOLNIR_FORCEFIELD_GLOBAL_PAIR_INTEARACTION_HPP
#define MJOLNIR_FORCEFIELD_GLOBAL_PAIR_INTEARACTION_HPP
#include <mjolnir/forcefield/global/ParameterList.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/GlobalInteractionBase.hpp>
#include <mjolnir/core/SpatialPartitionBase.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/string.hpp>
#include <memory>

namespace mjolnir
{

// Normal pair interaction.
//
// For the following potentials, there are specialized definitions to speedup
// calculation.
// - ExcludedVolume
//   - GlobalPairExcludedVolumeInteraction.hpp
// - LennardJones
//   - GlobalPairLennardJonesInteraction.hpp
// By combining potential calculations, we can omit `sqrt()` that is usually
// used to calculate the distance between particles.
//
template<typename traitsT, typename potentialT>
class GlobalPairInteraction final : public GlobalInteractionBase<traitsT>
{
  public:
    using traits_type         = traitsT;
    using potential_type      = potentialT;
    using base_type           = GlobalInteractionBase<traitsT>;
    using real_type           = typename base_type::real_type;
    using coordinate_type     = typename base_type::coordinate_type;
    using system_type         = typename base_type::system_type;
    using topology_type       = typename base_type::topology_type;
    using boundary_type       = typename base_type::boundary_type;
    using partition_type      = SpatialPartition<traits_type, potential_type>;
    using parameter_list_type = ParameterList<traits_type, potential_type>;

  public:

    GlobalPairInteraction(potential_type&& pot, parameter_list_type&& para,
                          partition_type&& part)
        : potential_(std::move(pot)), parameters_(std::move(para)),
          partition_(std::move(part))
    {}
    ~GlobalPairInteraction() override {}

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());

        this->potential_ .initialize(sys);
        this->parameters_.initialize(sys, topol, potential_);
        this->partition_ .initialize(sys, parameters_.cref());
    }

    /*! @brief update parameters (e.g. temperature, ionic strength, ...)  *
     *  @details A method that change system parameters (e.g. Annealing), *
     *           the method is bound to call this function after changing *
     *           parameters.                                              */
    void update(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());

        this->potential_ .update(sys);
        this->parameters_.update(sys, topol, potential_);
        this->partition_ .initialize(sys, this->parameters_.cref());
    }

    void reduce_margin(const real_type dmargin, const system_type& sys) override
    {
        this->partition_.reduce_margin(dmargin, sys, this->parameters_.cref());
        return;
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        this->partition_.scale_margin(scale, sys, this->parameters_.cref());
        return;
    }

    void calc_force (system_type& sys)        const noexcept override
    {
        this->template calc_force_virial_energy_impl<false, false>(sys);
        return;
    }
    void calc_force_and_virial(system_type& sys) const noexcept override
    {
        this->template calc_force_virial_energy_impl<false, true>(sys);
        return;
    }
    real_type calc_force_and_energy(system_type& sys) const noexcept override
    {
        return this->template calc_force_virial_energy_impl<true, false>(sys);
    }
    real_type calc_force_virial_energy(system_type& sys) const noexcept override
    {
        return this->template calc_force_virial_energy_impl<true, true>(sys);
    }

    real_type calc_energy(const system_type&)     const noexcept override;

    std::string name() const override
    {
        return "Pair:"_s + potential_type::name();
    }

    parameter_list_type const& parameters() const noexcept {return parameters_.cref();}
    parameter_list_type &      parameters()       noexcept {return parameters_.cref();}

    potential_type const& potential() const noexcept {return potential_;}
    potential_type&       potential()       noexcept {return potential_;}

    partition_type const& partition() const noexcept {return partition_;}

    // deep-copy through the base class, GlobalInteractionBase*
    base_type* clone() const override
    {
        return new GlobalPairInteraction(*this);
    }

  private:

    template<bool NeedEnergy, bool NeedVirial>
    real_type calc_force_virial_energy_impl(system_type& sys) const noexcept;

  private:

    potential_type      potential_;
    parameter_list_type parameters_;
    partition_type      partition_;

#ifdef MJOLNIR_WITH_OPENMP
    // OpenMP implementation uses its own implementation to run it in parallel.
    // So this implementation should not be instanciated with OpenMP Traits.
    static_assert(!is_openmp_simulator_traits<traits_type>::value,
                  "this is the default implementation, not for OpenMP");
#endif
};

template<typename traitsT, typename potT>
template<bool NeedEnergy, bool NeedVirial>
GlobalPairInteraction<traitsT, potT>::real_type
GlobalPairInteraction<traitsT, potT>::calc_force_virial_energy_impl(
        system_type& sys) const noexcept
{
    real_type energy = 0.0;
    const auto leading_participants = this->parameters_.leading_participants();
    for(std::size_t idx=0; idx<leading_participants.size(); ++idx)
    {
        const auto i = leading_participants[idx];
        for(const auto& ptnr : this->partition_.partners(i))
        {
            const auto  j    = ptnr.index;
            const auto& para = ptnr.parameter();

            const auto rij =
                sys.adjust_direction(sys.position(i), sys.position(j));
            const real_type l2 = math::length_sq(rij); // |rij|^2
            const real_type rl = math::rsqrt(l2);      // 1 / |rij|
            const real_type l  = l2 * rl;              // |rij|^2 / |rij|
            const real_type f_mag = potential_.derivative(l, para);

            // if length exceeds cutoff, potential returns just 0.
            if(f_mag == 0.0){continue;}

            const coordinate_type f = rij * (f_mag * rl);
            sys.force(i) += f;
            sys.force(j) -= f;

            if(NeedVirial) // (rj - ri) * Fj = (ri - rj) * Fi
            {
                sys.virial() += math::tensor_product(rij, -f);
            }
            if(NeedEnergy)
            {
                energy += potential_.potential(l, para);
            }
        }
    }
    return energy;
}

template<typename traitsT, typename potT>
typename GlobalPairInteraction<traitsT, potT>::real_type
GlobalPairInteraction<traitsT, potT>::calc_energy(
        const system_type& sys) const noexcept
{
    real_type E = 0.0;

    const auto leading_participants = this->parameters_.leading_participants();
    for(std::size_t idx=0; idx<leading_participants.size(); ++idx)
    {
        const auto i = leading_participants[idx];
        for(const auto& ptnr : this->partition_.partners(i))
        {
            const auto  j    = ptnr.index;
            const auto& para = ptnr.parameter();

            const real_type l = math::length(
                sys.adjust_direction(sys.position(i), sys.position(j)));
            E += potential_.potential(l, para);
        }
    }
    return E;
}

} // mjolnir

#include <mjolnir/forcefield/global/GlobalPairExcludedVolumeInteraction.hpp>
#include <mjolnir/forcefield/global/GlobalPairLennardJonesInteraction.hpp>

// #ifdef MJOLNIR_SEPARATE_BUILD
// // explicitly specialize BondAngleInteraction with LocalPotentials
// #include <mjolnir/core/BoundaryCondition.hpp>
// #include <mjolnir/core/SimulatorTraits.hpp>
// #include <mjolnir/forcefield/global/DebyeHuckelPotential.hpp>
// 
// namespace mjolnir
// {
// 
// // EXV, L-J and UL-J have its own specialization, so DO NOT specialize here.
// 
// // D-H
// extern template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        DebyeHuckelPotential<SimulatorTraits<double, UnlimitedBoundary>       >>;
// extern template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        DebyeHuckelPotential<SimulatorTraits<float,  UnlimitedBoundary>       >>;
// extern template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
// extern template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
// } // mjolnir
// #endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_GLOBAL_PAIR_INTEARACTION */
