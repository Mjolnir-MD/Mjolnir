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
    using parameter_type      = real_type;
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
    GlobalStoichiometricInteractionPotential(
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {}
    ~GlobalStoichiometricInteractionPotential() = default;
    GlobalStoichiometricInteractionPotential(const GlobalStoichiometricInteractionPotential&) = default;
    GlobalStoichiometricInteractionPotential(GlobalStoichiometricInteractionPotential&&) = default;
    GlobalStoichiometricInteractionPotential& operator=(const GlobalStoichiometricInteractionPotential&) = default;
    GlobalStoichiometricInteractionPotential& operator=(GlobalStoichiometricInteractionPotential&&) = default;

    // this value will be stored in NeighborList.
    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept
    {
        return pair_parameter_type(0);
    }

    real_type max_cutoff_length() const
    {
        return *(std::max_element(parameters_.cbegin(), parameters_.cend()));
    }

    // -----------------------------------------------------------------------
    // for spatial partitions
    //
    // Here, the default implementation uses Newton's 3rd law to reduce
    // calculation. For an interacting pair (i, j), forces applied to i and j
    // are equal in magnitude and opposite in direction. So, if a pair (i, j) is
    // listed, (j, i) is not needed.
    //     See implementation of VerletList, CellList and GlobalPairInteraction
    // for more details about the usage of these functions.

    std::vector<std::size_t> const& participants() const noexcept {return participants_;}

    range<typename std::vector<std::size_t>::const_iterator>
    leading_participants() const noexcept
    {
        return make_range(participants_.begin(), std::prev(participants_.end()));
    }
    range<typename std::vector<std::size_t>::const_iterator>
    possible_partners_of(const std::size_t participant_idx,
                         const std::size_t /*particle_idx*/) const noexcept
    {
        return make_range(participants_.begin() + participant_idx + 1,
                          participants_.end());
    }
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept
    {
        return (i < j) && !exclusion_list_.is_excluded(i, j);
    }

    // used by Observer.
    static const char* name() noexcept {return "ExcludedVolume";}

  private:

    std::vector<parameter_type> parameters_;
    std::vector<std::size_t>         participants_;

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
