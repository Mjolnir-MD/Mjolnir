#ifndef MJOLNIR_FORCEFIELD_GLOBAL_STOICHIOMETRIC_INTERACTION_POTENTIAL_HPP
#define MJOLNIR_FORCEFIELD_GLOBAL_STOICHIOMETRIC_INTERACTION_POTENTIAL_HPP
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>

namespace mjolnir
{

template<typename traitsT>
class GlobalStoichiometricInteractionPotential
{
  public:

    using traits_type         = traitsT;
    using real_type           = typename traits_type::real_type;
    using system_type         = System<traits_type>;
    using parameter_type      = real_type; // \delta r
    using pair_parameter_type = parameter_type;

    // ------------------------------------------------------------------------
    // topology stuff
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList <traits_type>;

  public:
    GlobalStoichiometricInteractionPotential(real_type v0,
        std::vector<std::size_t>&& first_kind_participants,
        std::vector<std::size_t>&& second_kind_participants,
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : v0_(v0), inv_v0_(1.0 / v0),
          participants_(first_kind_participants),
          first_kind_participants_num_(first_kind_participants.size()),
          second_kind_participants_num_(second_kind_participants.size()),
          exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {
        MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
        MJOLNIR_LOG_FUNCTION_DEBUG();
        participants_.insert(participants_.end(),
                             second_kind_participants.begin(), second_kind_participants.end());
    }
    ~GlobalStoichiometricInteractionPotential() = default;
    GlobalStoichiometricInteractionPotential(const GlobalStoichiometricInteractionPotential&) = default;
    GlobalStoichiometricInteractionPotential(GlobalStoichiometricInteractionPotential&&) = default;
    GlobalStoichiometricInteractionPotential& operator=(const GlobalStoichiometricInteractionPotential&) = default;
    GlobalStoichiometricInteractionPotential& operator=(GlobalStoichiometricInteractionPotential&&) = default;

    pair_parameter_type prepare_params(std::size_t, std::size_t) const noexcept
    {
        return v0_;
    }

    real_type potential(const real_type r) const noexcept
    {
        if(v0_ < r){return 0.0;}

        const real_type r_inv_v0   = r * inv_v0_;
        const real_type r_inv_v0_2 = r_inv_v0 * r_inv_v0;
        const real_type r_inv_v0_3 = r_inv_v0_2 * r_inv_v0;
        return 1.0 - 3.0 * r_inv_v0_2 + 2.0 * r_inv_v0_3;
    }

    real_type derivative(const real_type r) const noexcept
    {
        if(v0_ < r){return 0.0;}

        const real_type r_inv_v0   = r * inv_v0_;
        const real_type r_inv_v0_2 = r_inv_v0 * r_inv_v0;
        return 6.0 * (r_inv_v0_2 - r_inv_v0);
    }

    std::size_t first_kind_participants_num()  const noexcept {return first_kind_participants_num_;}
    std::size_t second_kind_participants_num() const noexcept {return second_kind_participants_num_;}

    real_type max_cutoff_length() const noexcept
    {
        return 2.0 * v0_;
    }

    void initialize(const system_type& sys, const topology_type& topol) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->update(sys, topol);
        return;
    }

    void update(const system_type& sys, const topology_type& topol) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // update exclusion list based on sys.topology()
        exclusion_list_.make(sys, topol);
    }

    // -----------------------------------------------------------------------
    // for spatial partitions

    const std::vector<std::size_t>& participants() const noexcept {return participants_;}

    range<typename std::vector<std::size_t>::const_iterator>
    leading_participants() const noexcept
    {
        return make_range(participants_.begin(), 
                          participants_.begin() + first_kind_participants_num_);
    }
    range<typename std::vector<std::size_t>::const_iterator>
    possible_partners_of(const std::size_t /*participant_idx*/,
                         const std::size_t /*particle_idx*/) const noexcept
    {
        return make_range(participants_.begin() + first_kind_participants_num_,
                          participants_.end());
    }
    range<typename std::vector<std::size_t>::const_iterator>
    following_participants() const noexcept
    {
        return possible_partners_of(/*dummy param*/0, /*dummy param*/0);
    }
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept
    {
        return !exclusion_list_.is_excluded(i, j);
    }

    // -----------------------------------------------------------------------
    // used by Observer.
    static const char* name() noexcept {return "GlobalStoichiometricInteraction";}

  private:

    real_type                v0_;
    real_type                inv_v0_;
    std::vector<std::size_t> participants_;
    std::size_t              first_kind_participants_num_;
    std::size_t              second_kind_participants_num_;

    exclusion_list_type exclusion_list_;
};

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize major use-cases
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

namespace mjolnir
{
extern template class GlobalStoichiometricInteractionPotential<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class GlobalStoichiometricInteractionPotential<SimulatorTraits<float , UnlimitedBoundary>       >;
extern template class GlobalStoichiometricInteractionPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class GlobalStoichiometricInteractionPotential<SimulatorTraits<float , CuboidalPeriodicBoundary>>;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif // MJOLNIR_POTENTIAL_GLOBAL_STOICHIOMETRIC_INTERACTION_POTENTIAL_HPP
