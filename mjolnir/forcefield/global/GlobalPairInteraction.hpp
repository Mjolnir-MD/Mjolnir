#ifndef MJOLNIR_INTEARACTION_GLOBAL_PAIR_INTEARACTION_HPP
#define MJOLNIR_INTEARACTION_GLOBAL_PAIR_INTEARACTION_HPP
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
// - UniformLennardJones
//   - GlobalPairUniformLennardJonesInteraction.hpp
// By combining potential calculations, we can omit `sqrt()` that is usually
// used to calculate the distance between particles.
//
template<typename traitsT, typename potentialT>
class GlobalPairInteraction final : public GlobalInteractionBase<traitsT>
{
  public:
    using traits_type     = traitsT;
    using potential_type  = potentialT;
    using base_type       = GlobalInteractionBase<traitsT>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using topology_type   = typename base_type::topology_type;
    using boundary_type   = typename base_type::boundary_type;
    using partition_type  = SpatialPartition<traits_type, potential_type>;

  public:
    GlobalPairInteraction()           = default;
    ~GlobalPairInteraction() override {}

    GlobalPairInteraction(potential_type&& pot, partition_type&& part)
        : potential_(std::move(pot)), partition_(std::move(part))
    {}

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());
        this->potential_.initialize(sys, topol);
        this->partition_.initialize(sys, this->potential_);
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
        this->potential_.update(sys, topol);
        // potential update may change the cutoff length!
        this->partition_.initialize(sys, this->potential_);
    }

    void reduce_margin(const real_type dmargin, const system_type& sys) override
    {
        this->partition_.reduce_margin(dmargin, sys, this->potential_);
        return;
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        this->partition_.scale_margin(scale, sys, this->potential_);
        return;
    }

    void      calc_force (system_type&)           const noexcept override;
    real_type calc_energy(const system_type&)     const noexcept override;
    real_type calc_force_and_energy(system_type&) const noexcept override;

    std::string name() const override
    {return "Pair:"_s + potential_type::name();}

    potential_type const& potential() const noexcept {return potential_;}
    potential_type &      potential()       noexcept {return potential_;}

    base_type* clone() const override
    {
        return new GlobalPairInteraction(
                potential_type(potential_), partition_type(partition_));
    }

  private:

    potential_type potential_;
    partition_type partition_;

#ifdef MJOLNIR_WITH_OPENMP
    // OpenMP implementation uses its own implementation to run it in parallel.
    // So this implementation should not be instanciated with OpenMP Traits.
    static_assert(!is_openmp_simulator_traits<traits_type>::value,
                  "this is the default implementation, not for OpenMP");
#endif
};

template<typename traitsT, typename potT>
void GlobalPairInteraction<traitsT, potT>::calc_force(
        system_type& sys) const noexcept
{
    const auto leading_participants = this->potential_.leading_participants();
    for(std::size_t idx=0; idx<leading_participants.size(); ++idx)
    {
        const auto i = leading_participants[idx];
        for(const auto& ptnr : this->partition_.partners(i))
        {
            const auto  j     = ptnr.index;
            const auto& param = ptnr.parameter();

            const auto rij =
                sys.adjust_direction(sys.position(i), sys.position(j));
            const real_type l2 = math::length_sq(rij); // |rij|^2
            const real_type rl = math::rsqrt(l2);      // 1 / |rij|
            const real_type l  = l2 * rl;              // |rij|^2 / |rij|
            const real_type f_mag = potential_.derivative(l, param);

            // if length exceeds cutoff, potential returns just 0.
            if(f_mag == 0.0){continue;}

            const coordinate_type f = rij * (f_mag * rl);
            sys.force(i) += f;
            sys.force(j) -= f;
        }
    }
    return ;
}

template<typename traitsT, typename potT>
typename GlobalPairInteraction<traitsT, potT>::real_type
GlobalPairInteraction<traitsT, potT>::calc_energy(
        const system_type& sys) const noexcept
{
    real_type E = 0.0;

    const auto leading_participants = this->potential_.leading_participants();
    for(std::size_t idx=0; idx<leading_participants.size(); ++idx)
    {
        const auto i = leading_participants[idx];
        for(const auto& ptnr : this->partition_.partners(i))
        {
            const auto  j     = ptnr.index;
            const auto& param = ptnr.parameter();

            const real_type l = math::length(
                sys.adjust_direction(sys.position(i), sys.position(j)));
            E += potential_.potential(l, param);
        }
    }
    return E;
}

template<typename traitsT, typename potT>
typename GlobalPairInteraction<traitsT, potT>::real_type
GlobalPairInteraction<traitsT, potT>::calc_force_and_energy(
        system_type& sys) const noexcept
{
    real_type energy = 0.0;
    const auto leading_participants = this->potential_.leading_participants();
    for(std::size_t idx=0; idx<leading_participants.size(); ++idx)
    {
        const auto i = leading_participants[idx];
        for(const auto& ptnr : this->partition_.partners(i))
        {
            const auto  j     = ptnr.index;
            const auto& param = ptnr.parameter();

            const auto rij =
                sys.adjust_direction(sys.position(i), sys.position(j));
            const real_type l2 = math::length_sq(rij); // |rij|^2
            const real_type rl = math::rsqrt(l2);      // 1 / |rij|
            const real_type l  = l2 * rl;              // |rij|^2 / |rij|
            const real_type f_mag = potential_.derivative(l, param);

            // if length exceeds cutoff, potential returns just 0.
            if(f_mag == 0.0){continue;}

            energy += potential_.potential(l, param);
            const coordinate_type f = rij * (f_mag * rl);
            sys.force(i) += f;
            sys.force(j) -= f;
        }
    }
    return energy;
}


} // mjolnir

#include <mjolnir/forcefield/global/GlobalPairExcludedVolumeInteraction.hpp>
#include <mjolnir/forcefield/global/GlobalPairLennardJonesInteraction.hpp>
#include <mjolnir/forcefield/global/GlobalPairUniformLennardJonesInteraction.hpp>

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize BondAngleInteraction with LocalPotentials
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/global/DebyeHuckelPotential.hpp>

namespace mjolnir
{

// EXV, L-J and UL-J have its own specialization, so DO NOT specialize here.

// D-H
extern template class GlobalPairInteraction<SimulatorTraits<double, UnlimitedBoundary>,        DebyeHuckelPotential<SimulatorTraits<double, UnlimitedBoundary>       >>;
extern template class GlobalPairInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        DebyeHuckelPotential<SimulatorTraits<float,  UnlimitedBoundary>       >>;
extern template class GlobalPairInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class GlobalPairInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_GLOBAL_PAIR_INTEARACTION */
