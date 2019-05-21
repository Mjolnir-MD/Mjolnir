#ifndef MJOLNIR_OMP_GLOBAL_PAIR_LENNARD_JONES_INTEARACTION_HPP
#define MJOLNIR_OMP_GLOBAL_PAIR_LENNARD_JONES_INTEARACTION_HPP
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/interaction/global/GlobalPairInteraction.hpp>
#include <mjolnir/potential/global/DebyeHuckelPotential.hpp>

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT,
         typename partitionT>
class GlobalPairInteraction<
    OpenMPSimulatorTraits<realT, boundaryT>,
    DebyeHuckelPotential<realT>,
    partitionT
    > final : public GlobalInteractionBase<OpenMPSimulatorTraits<realT, boundaryT>>
{
  public:

    using traits_type     = OpenMPSimulatorTraits<realT, boundaryT>;
    using potential_type  = DebyeHuckelPotential<realT>;
    using partition_type  = partitionT;
    using base_type       = GlobalInteractionBase<traits_type>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using boundary_type   = typename base_type::boundary_type;

  public:
    GlobalPairInteraction()  = default;
    ~GlobalPairInteraction() = default;

    GlobalPairInteraction(potential_type&& pot, partition_type&& part)
        : potential_(std::move(pot)), partition_(std::move(part))
    {}

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());
        this->potential_.initialize(sys);
        this->partition_.initialize(sys, this->potential_);
    }

    /*! @brief update parameters (e.g. temperature, ionic strength, ...)  *
     *  @details A method that change system parameters (e.g. Annealing), *
     *           the method is bound to call this function after changing *
     *           parameters.                                              */
    void update(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());
        this->potential_.update(sys);
        // potential update may change the cutoff length!
        this->partition_.initialize(sys, this->potential_);
    }

    void update_margin(const real_type dmargin, const system_type& sys) override
    {
        this->partition_.update(dmargin, sys, this->potential_);
        return;
    }

    void calc_force (system_type& sys) const noexcept override
    {
#pragma omp for nowait
        for(std::size_t idx=0; idx < this->potential_.participants().size(); ++idx)
        {
            const auto i = this->potential_.participants()[idx];
            for(const auto& ptnr : this->partition_.partners(i))
            {
                const auto  j     = ptnr.index;
                const auto& param = ptnr.parameter();

                const coordinate_type rij =
                    sys.adjust_direction(sys.position(j) - sys.position(i));
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
#pragma omp parallel for reduction(+:E)
        for(std::size_t idx=0; idx < this->potential_.participants().size(); ++idx)
        {
            const auto i = this->potential_.participants()[idx];
            for(const auto& ptnr : this->partition_.partners(i))
            {
                const auto  j     = ptnr.index;
                const auto& param = ptnr.parameter();

                const real_type l = math::length(
                    sys.adjust_direction(sys.position(j) - sys.position(i)));
                E += potential_.potential(l, param);
            }
        }
        return E;
    }

    std::string name() const override
    {return "Pair:"_s + potential_type::name();}

  private:

    potential_type potential_;
    partition_type partition_;
};

} // mjolnir
#endif /* MJOLNIR_GLOBAL_PAIR_INTEARACTION */
