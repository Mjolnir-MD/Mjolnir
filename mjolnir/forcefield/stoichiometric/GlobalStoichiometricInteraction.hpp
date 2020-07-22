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

    using traits_type     = traitsT;
    using base_type       = GlobalInteractionBase<traits_type>;
    using index_type      = std::size_t;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using topology_type   = typename base_type::topology_type;
    using potential_type  = GlobalStoichiometricInteractionPotential<traits_type>;
    using partition_type  = SpatialPartition<traitsT, potential_type>;

  public:
    GlobalStoichiometricInteraction(potential_type&& pot, partition_type&& part,
        real_type epsilon,
        std::size_t first_coef, std::size_t second_coef)
        : potential_(std::move(pot)), partition_(std::move(part)),
          epsilon_(epsilon),
          first_coef_(first_coef),
          second_coef_(second_coef)
    {
        std::size_t first_kind_particles_num  = potential_.first_kind_participants_num();
        std::size_t second_kind_particles_num = potential_.second_kind_participants_num();
        pot_derivs_buff_.resize(first_kind_particles_num);
        potentials_buff_.resize(first_kind_particles_num);
        for(std::size_t idx = 0; idx<first_kind_particles_num; ++idx)
        {
            potentials_buff_[idx].resize(second_kind_particles_num);
            pot_derivs_buff_[idx].resize(second_kind_particles_num);
        }
        pot_sum_for_first_       .resize(first_kind_particles_num );
        pot_sum_for_second_      .resize(second_kind_particles_num);
        pot_deriv_sum_for_first_ .resize(first_kind_particles_num );
        pot_deriv_sum_for_second_.resize(second_kind_particles_num);
        activated_for_first_     .resize(first_kind_particles_num );
        activated_for_second_    .resize(second_kind_particles_num);
        act_deriv_for_first_     .resize(first_kind_particles_num );
        act_deriv_for_second_    .resize(second_kind_particles_num);
    }
    ~GlobalStoichiometricInteraction() {}

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());
        potential_.initialize(sys, topol);
        partition_.initialize(sys, this->potential_);
    }

    void update(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());
        this->potential_.update(sys,topol);
        this->partition_.initialize(sys, potential_);
    }

    void reduce_margin(const real_type dmargin, const system_type& sys) override
    {
        partition_.reduce_margin(dmargin, sys, potential_);
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        partition_.scale_margin(scale, sys, potential_);
    }

    void      calc_force (system_type&)       const noexcept override;
    real_type calc_energy(const system_type&) const noexcept override;

    std::string name() const override {return "Stoichiometric";}

    potential_type const& potential() const noexcept {return potential_;}
    potential_type &      potential()       noexcept {return potential_;}

    base_type* clone() const override
    {
        return new GlobalStoichiometricInteraction(
                potential_type(potential_), partition_type(partition_),
                epsilon_, first_coef_, second_coef_);
    }

  private:
    real_type activation_func(const real_type x) const noexcept
    {
        // This function correspond to
        //                     0             (x < -1)
        // sigma(x) =  1/4 (2 + 3x - x^3)    (-1 <= x <= 1)
        //                     1             (1 < x)

        if(x < -1.0){return 0;}
        if(1.0 < x) {return 1;}
        const real_type x_2 = x * x;
        const real_type x_3 = x_2 * x;
        return (2.0 + 3.0 * x - x_3) * 0.25;
    }

    real_type deriv_activation_func(const real_type x) const
    {
        if(x < -1.0 || 1.0 < x){return 0;}
        return 0.75 - 0.75 * x * x;
    }

  private:

    potential_type                      potential_;
    partition_type                      partition_;

    real_type                           epsilon_;
    std::size_t                         first_coef_;
    std::size_t                         second_coef_;


    // -----------------------------------------------------------------------
    // Variables for buffering intermediate value.
    // These value can change in calc_force and calc_energy function.
    mutable std::vector<std::vector<real_type>>       potentials_buff_;
    mutable std::vector<std::vector<coordinate_type>> pot_derivs_buff_;
    mutable std::vector<real_type>                    pot_sum_for_first_;
    mutable std::vector<real_type>                    pot_sum_for_second_;
    // sum of derivation of potential function for specific first particle.
    mutable std::vector<coordinate_type>              pot_deriv_sum_for_first_;
    // sum of derivation of potential function for specific second particle.
    mutable std::vector<coordinate_type>              pot_deriv_sum_for_second_;
    // buffer container for pre calculated value of activation function for each particle.
    mutable std::vector<real_type>                    activated_for_first_;
    mutable std::vector<real_type>                    activated_for_second_;
    // buffer container for pre calculated derivation value of activation function
    // for each particle.
    mutable std::vector<real_type>                    act_deriv_for_first_;
    mutable std::vector<real_type>                    act_deriv_for_second_;
    // temporary force container for calc_force.
    mutable coordinate_type                           tmp_coord_first_;
    mutable coordinate_type                           tmp_coord_second_;
};

template<typename traitsT>
void GlobalStoichiometricInteraction<traitsT>::calc_force(system_type& sys) const noexcept
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();

    std::fill(pot_sum_for_first_ .begin(), pot_sum_for_first_ .end(), 0.0);
    std::fill(pot_sum_for_second_.begin(), pot_sum_for_second_.end(), 0.0);
    std::fill(pot_deriv_sum_for_first_ .begin(), pot_deriv_sum_for_first_ .end(),
              math::make_coordinate<coordinate_type>(0.0, 0.0, 0.0));
    std::fill(pot_deriv_sum_for_second_.begin(), pot_deriv_sum_for_second_.end(),
              math::make_coordinate<coordinate_type>(0.0, 0.0, 0.0));
    const std::size_t first_kind_participants_num  = potential_.first_kind_participants_num();
    const std::size_t second_kind_participants_num = potential_.second_kind_participants_num();

    const auto leading_participants   = potential_.leading_participants();
    const auto following_participants = potential_.following_participants();
    // pre calculation for pot_sum and pot_deriv_sum for each particle.
    for(std::size_t idx_first=0; idx_first<first_kind_participants_num; ++idx_first)
    {
        const index_type i = leading_participants[idx_first];
        std::vector<coordinate_type>& derivs_buff_first = pot_derivs_buff_[idx_first];
        std::vector<real_type>&       pots_buff_first   = potentials_buff_[idx_first];
        for(std::size_t idx_second=0; idx_second<second_kind_participants_num; ++idx_second)
        {
            const index_type j = following_participants[idx_second];
            const coordinate_type rij =
               sys.adjust_direction(sys.position(j) - sys.position(i));
            const real_type       l2    = math::length_sq(rij); // |rij|^2
            const real_type       rl    = math::rsqrt(l2);      // 1 / |rij|
            const coordinate_type e     = rij * rl;             // rij / |rij|
            const real_type       l     = l2 * rl;
            const coordinate_type deriv = potential_.derivative(l) * e;
            const real_type       pot   = potential_.potential(l);
            derivs_buff_first[idx_second]         =  deriv;
            pot_deriv_sum_for_first_[idx_first]   += deriv;
            pot_deriv_sum_for_second_[idx_second] -= deriv;
            pots_buff_first[idx_second]     =  pot;
            pot_sum_for_first_[idx_first]   += pot;
            pot_sum_for_second_[idx_second] += pot;
        }
    }

    // pre calculation for the value of activated function for each particle.
    for(std::size_t idx_first=0; idx_first<first_kind_participants_num; ++idx_first)
    {
        const real_type  x = 2.0 * (first_coef_ + 0.5 - pot_sum_for_first_[idx_first]);
        activated_for_first_[idx_first] = activation_func(x);
        act_deriv_for_first_[idx_first] = deriv_activation_func(x);
    }
    for(std::size_t idx_second=0; idx_second<second_kind_participants_num; ++idx_second)
    {
        const real_type  x = 2.0 * (second_coef_ + 0.5 - pot_sum_for_second_[idx_second]);
        activated_for_second_[idx_second] = activation_func(x);
        act_deriv_for_second_[idx_second] = deriv_activation_func(x);
    }

    // force calculation
    for(std::size_t idx_first=0; idx_first<first_kind_participants_num; ++idx_first)
    {
        const std::vector<real_type>&       pots_buff_first       = potentials_buff_[idx_first];
        const std::vector<coordinate_type>& pot_derivs_buff_first = pot_derivs_buff_[idx_first];
        for(std::size_t idx_second=0; idx_second<second_kind_participants_num; ++idx_second)
        {
            const real_type pots_buff_first_second = pots_buff_first[idx_second];
            const real_type activated_for_first    = activated_for_first_[idx_first];
            const real_type activated_for_second   = activated_for_second_[idx_second];
            const coordinate_type& pot_derivs_buff_first_second = pot_derivs_buff_first[idx_second];
            const real_type term1_coef =
                activated_for_second * pots_buff_first_second * act_deriv_for_first_[idx_first];
            const real_type term2_coef =
                activated_for_first  * pots_buff_first_second * act_deriv_for_second_[idx_second];
            const real_type term3_coef = activated_for_first * activated_for_second;

            const index_type i = leading_participants[idx_first];
            const index_type j = following_participants[idx_second];
            sys.force(i) += epsilon_ *
                            (2.0 * term1_coef * pot_deriv_sum_for_first_[idx_first] +
                            (2.0 * term2_coef - term3_coef) * pot_derivs_buff_first_second);
            sys.force(j) += epsilon_ *
                            (2.0 * term2_coef * pot_deriv_sum_for_second_[idx_second] +
                            (2.0 * term1_coef - term3_coef) * -pot_derivs_buff_first_second);
        }
    }

    return;
}

template<typename traitsT>
typename GlobalStoichiometricInteraction<traitsT>::real_type
GlobalStoichiometricInteraction<traitsT>::calc_energy(const system_type& sys) const noexcept
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();

    const std::size_t first_kind_participants_num  = potential_.first_kind_participants_num();
    const std::size_t second_kind_participants_num = potential_.second_kind_participants_num();
    const auto leading_participants   = potential_.leading_participants();
    const auto following_participants = potential_.following_participants();

    std::fill(pot_sum_for_first_ .begin(), pot_sum_for_first_ .end(), 0.0);
    std::fill(pot_sum_for_second_.begin(), pot_sum_for_second_.end(), 0.0);
    for(std::size_t idx_first=0; idx_first<first_kind_participants_num; ++idx_first)
    {
        const index_type i = leading_participants[idx_first];
        std::vector<real_type>& pots_buff_first   = potentials_buff_  [idx_first];
        real_type&              pot_sum_for_first = pot_sum_for_first_[idx_first];
        for(std::size_t idx_second=0; idx_second<second_kind_participants_num; ++idx_second)
        {
            const index_type j = following_participants[idx_second];
            const coordinate_type rij =
                sys.adjust_direction(sys.position(j) - sys.position(i));
            const real_type l2  = math::length_sq(rij); // |rij|^2
            const real_type rl  = math::rsqrt(l2);      // 1 / |rij|
            const real_type l   = l2 * rl;              // |rij|
            const real_type pot = potential_.potential(l);
            pots_buff_first    [idx_second] =  pot;
            pot_sum_for_first               += pot;
            pot_sum_for_second_[idx_second] += pot;
        }
    }

    real_type retval(0.0);
    for(std::size_t idx_first=0; idx_first<first_kind_participants_num; ++idx_first) {
        const std::vector<real_type>& pots_buff_first = potentials_buff_[idx_first];
        for(std::size_t idx_second=0; idx_second<second_kind_participants_num; ++idx_second)
        {
            const real_type x_for_first  =
                2.0 * (first_coef_  + 0.5 - pot_sum_for_first_[idx_first]);
            const real_type x_for_second =
                2.0 * (second_coef_ + 0.5 - pot_sum_for_second_[idx_second]);
            retval += activation_func(x_for_first) * activation_func(x_for_second) *
                      pots_buff_first[idx_second];
        }
    }
    retval *= -epsilon_;

    return retval;
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
