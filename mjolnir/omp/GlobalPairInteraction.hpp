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
    using parameter_list_type = ParameterList<traits_type, potential_type>;

  public:
    GlobalPairInteraction()  = default;
    ~GlobalPairInteraction() override {}

    GlobalPairInteraction(potential_type&& pot, parameter_list_type&& para,
                          partition_type&& part)
        : potential_(std::move(pot)), parameters_(std::move(para)),
          partition_(std::move(part))
    {}

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

    void calc_force (system_type& sys) const noexcept override
    {
        this->template calc_force_and_virial_impl<false>(sys);
        return ;
    }
    void calc_force_and_virial(system_type& sys) const noexcept override
    {
        this->template calc_force_and_virial_impl<true>(sys);
        return ;
    }

    real_type calc_force_and_energy(system_type& sys) const noexcept override
    {
        return this->template calc_force_virial_energy_impl<false>(sys);
    }
    real_type calc_force_virial_energy(system_type& sys) const noexcept override
    {
        return this->template calc_force_virial_energy_impl<true>(sys);
    }

    real_type calc_energy(const system_type& sys) const noexcept override
    {
        real_type E = 0.0;
        const auto leading_participants = this->parameters_.leading_participants();
#pragma omp parallel for reduction(+:E)
        for(std::size_t idx=0; idx < leading_participants.size(); ++idx)
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

    std::string name() const override
    {return "Pair:"_s + potential_type::name();}

    base_type* clone() const override
    {
        return new GlobalPairInteraction(*this);
    }

  private:

    template<bool NeedVirial>
    void calc_force_and_virial_impl(system_type& sys) const noexcept
    {
        const auto leading_participants = this->parameters_.leading_participants();
#pragma omp parallel for
        for(std::size_t idx=0; idx < leading_participants.size(); ++idx)
        {
            const std::size_t thread_id = omp_get_thread_num();
            const auto i = leading_participants[idx];
            for(const auto& ptnr : this->partition_.partners(i))
            {
                const auto  j    = ptnr.index;
                const auto& para = ptnr.parameter();

                const coordinate_type rij =
                    sys.adjust_direction(sys.position(i), sys.position(j));
                const real_type l2 = math::length_sq(rij); // |rij|^2
                const real_type rl = math::rsqrt(l2);      // 1 / |rij|
                const real_type l  = l2 * rl;              // |rij|^2 / |rij|
                const real_type f_mag = potential_.derivative(l, para);

                // if length exceeds cutoff, potential returns just 0.
                if(f_mag == 0.0){continue;}

                const coordinate_type f = rij * (f_mag * rl);
                sys.force_thread(thread_id, i) += f;
                sys.force_thread(thread_id, j) -= f;

                if(NeedVirial)
                {
                    sys.virial_thread(thread_id) += math::tensor_product(rij, -f);
                }
            }
        }
        return ;
    }

    template<bool NeedVirial>
    real_type calc_force_virial_energy_impl(system_type& sys) const noexcept
    {
        real_type energy = 0.0;
        const auto leading_participants = this->parameters_.leading_participants();
#pragma omp parallel for reduction(+:energy)
        for(std::size_t idx=0; idx < leading_participants.size(); ++idx)
        {
            const std::size_t thread_id = omp_get_thread_num();
            const auto i = leading_participants[idx];
            for(const auto& ptnr : this->partition_.partners(i))
            {
                const auto  j    = ptnr.index;
                const auto& para = ptnr.parameter();

                const coordinate_type rij =
                    sys.adjust_direction(sys.position(i), sys.position(j));
                const real_type l2 = math::length_sq(rij); // |rij|^2
                const real_type rl = math::rsqrt(l2);      // 1 / |rij|
                const real_type l  = l2 * rl;              // |rij|^2 / |rij|
                const real_type f_mag = potential_.derivative(l, para);

                // if length exceeds cutoff, potential returns just 0.
                if(f_mag == 0.0){continue;}

                energy += potential_.potential(l, para);

                const coordinate_type f = rij * (f_mag * rl);
                sys.force_thread(thread_id, i) += f;
                sys.force_thread(thread_id, j) -= f;

                if(NeedVirial)
                {
                    sys.virial_thread(thread_id) += math::tensor_product(rij, -f);
                }
            }
        }
        return energy;
    }


  private:

    potential_type      potential_;
    parameter_list_type parameters_;
    partition_type      partition_;
};

} // mjolnir

// specializations
#include <mjolnir/omp/GlobalPairExcludedVolumeInteraction.hpp>
#include <mjolnir/omp/GlobalPairLennardJonesInteraction.hpp>

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize BondAngleInteraction with LocalPotentials
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/forcefield/global/DebyeHuckelPotential.hpp>

namespace mjolnir
{

// EXV, L-J and UL-J have its own specialization, so DO NOT specialize here.

// D-H
extern template class GlobalPairInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>       , DebyeHuckelPotential<double>>;
extern template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       , DebyeHuckelPotential<float >>;
extern template class GlobalPairInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>>;
extern template class GlobalPairInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float >>;
} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_GLOBAL_PAIR_INTEARACTION */
