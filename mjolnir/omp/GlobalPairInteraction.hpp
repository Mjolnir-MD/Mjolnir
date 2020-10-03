#ifndef MJOLNIR_OMP_GLOBAL_PAIR_INTEARACTION_HPP
#define MJOLNIR_OMP_GLOBAL_PAIR_INTEARACTION_HPP
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/forcefield/global/GlobalPairInteraction.hpp>

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT,
         typename potentialT>
class GlobalPairInteraction<
    OpenMPSimulatorTraits<realT, boundaryT>, potentialT
    > final : public GlobalInteractionBase<OpenMPSimulatorTraits<realT, boundaryT>>
{
  public:
    using traits_type     = OpenMPSimulatorTraits<realT, boundaryT>;
    using potential_type  = potentialT;
    using base_type       = GlobalInteractionBase<traits_type>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using topology_type   = typename base_type::topology_type;
    using boundary_type   = typename base_type::boundary_type;
    using partition_type  = SpatialPartition<traits_type, potential_type>;

  public:
    GlobalPairInteraction()  = default;
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

    void calc_force (system_type& sys) const noexcept override
    {
        const auto leading_participants = this->potential_.leading_participants();
#pragma omp parallel for
        for(std::size_t idx=0; idx < leading_participants.size(); ++idx)
        {
            const auto i = leading_participants[idx];
            for(const auto& ptnr : this->partition_.partners(i))
            {
                const auto  j     = ptnr.index;
                const auto& param = ptnr.parameter();

                const coordinate_type rij =
                    sys.adjust_direction(sys.position(i), sys.position(j));
                const real_type l2 = math::length_sq(rij); // |rij|^2
                const real_type rl = math::rsqrt(l2);      // 1 / |rij|
                const real_type l  = l2 * rl;              // |rij|^2 / |rij|
                const real_type f_mag = potential_.derivative(l, param);

                // if length exceeds cutoff, potential returns just 0.
                if(f_mag == 0.0){continue;}

                const coordinate_type f = rij * (f_mag * rl);
                const std::size_t thread_id = omp_get_thread_num();
                sys.force_thread(thread_id, i) += f;
                sys.force_thread(thread_id, j) -= f;
            }
        }
        return ;
    }
    real_type calc_energy(const system_type& sys) const noexcept override
    {
        real_type E = 0.0;
        const auto leading_participants = this->potential_.leading_participants();
#pragma omp parallel for reduction(+:E)
        for(std::size_t idx=0; idx < leading_participants.size(); ++idx)
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

    real_type calc_force_and_energy(system_type& sys) const noexcept override
    {
        real_type energy = 0.0;
        const auto leading_participants = this->potential_.leading_participants();
#pragma omp parallel for reduction(+:energy)
        for(std::size_t idx=0; idx < leading_participants.size(); ++idx)
        {
            const auto i = leading_participants[idx];
            for(const auto& ptnr : this->partition_.partners(i))
            {
                const auto  j     = ptnr.index;
                const auto& param = ptnr.parameter();

                const coordinate_type rij =
                    sys.adjust_direction(sys.position(i), sys.position(j));
                const real_type l2 = math::length_sq(rij); // |rij|^2
                const real_type rl = math::rsqrt(l2);      // 1 / |rij|
                const real_type l  = l2 * rl;              // |rij|^2 / |rij|
                const real_type f_mag = potential_.derivative(l, param);

                // if length exceeds cutoff, potential returns just 0.
                if(f_mag == 0.0){continue;}

                energy += potential_.potential(l, param);

                const coordinate_type f = rij * (f_mag * rl);
                const std::size_t thread_id = omp_get_thread_num();
                sys.force_thread(thread_id, i) += f;
                sys.force_thread(thread_id, j) -= f;
            }
        }
        return energy;
    }

    std::string name() const override
    {return "Pair:"_s + potential_type::name();}

    base_type* clone() const override
    {
        return new GlobalPairInteraction(
                potential_type(potential_), partition_type(partition_));
    }


  private:

    potential_type potential_;
    partition_type partition_;
};

} // mjolnir

// specializations
#include <mjolnir/omp/GlobalPairExcludedVolumeInteraction.hpp>
#include <mjolnir/omp/GlobalPairLennardJonesInteraction.hpp>
#include <mjolnir/omp/GlobalPairUniformLennardJonesInteraction.hpp>

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize BondAngleInteraction with LocalPotentials
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/forcefield/global/DebyeHuckelPotential.hpp>

namespace mjolnir
{

// EXV, L-J and UL-J have its own specialization, so DO NOT specialize here.

// D-H
extern template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>       , DebyeHuckelPotential<OpenMPSimulatorTraits<double, UnlimitedBoundary>       >>;
extern template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       , DebyeHuckelPotential<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       >>;
extern template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_GLOBAL_PAIR_INTEARACTION */
