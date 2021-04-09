#ifndef MJOLNIR_FORCEFIELD_GLOBAL_STOICHIOMETRIC_INTERACTION_POTENTIAL_HPP
#define MJOLNIR_FORCEFIELD_GLOBAL_STOICHIOMETRIC_INTERACTION_POTENTIAL_HPP
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>

#include <algorithm>

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
    GlobalStoichiometricInteractionPotential(real_type v0, real_type range,
        std::vector<std::size_t>&& participants_a,
        std::vector<std::size_t>&& participants_b,
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : v0_(v0), v0_range_(v0 + range),
          inv_range_(1.0 / range), inv_range_6_(6.0 / range),
          participants_(participants_a),
          participants_a_(participants_a), participants_a_num_(participants_a.size()),
          participants_b_(participants_b),
          exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {
        MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
        MJOLNIR_LOG_FUNCTION_DEBUG();
        participants_.insert(participants_.end(),
                             participants_b.begin(), participants_b.end());
        std::sort(participants_.begin(), participants_.end());
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
        if(r < v0_){return 1.0;}
        else if(v0_range_ < r){return 0.0;}

        const real_type r_v0           = r - v0_;
        const real_type r_v0_inv_range = r_v0 * inv_range_;
        const real_type r_v0_inv_range_2  = r_v0_inv_range * r_v0_inv_range;
        const real_type r_v0_inv_range_3  = r_v0_inv_range_2 * r_v0_inv_range;
        return 1.0 - 3.0 * r_v0_inv_range_2 + 2.0 * r_v0_inv_range_3;
    }

    real_type derivative(const real_type r) const noexcept
    {
        if(r < v0_ || v0_range_ < r){return 0.0;}

        const real_type r_v0          = r - v0_;
        const real_type r_v0_inv_range   = r_v0 * inv_range_;
        const real_type r_v0_inv_range_2 = r_v0_inv_range * r_v0_inv_range;
        return inv_range_6_ * (r_v0_inv_range_2 - r_v0_inv_range);
    }

    std::vector<std::size_t> participants_a() const noexcept {return participants_a_;}
    std::vector<std::size_t> participants_b() const noexcept {return participants_b_;}

    std::size_t participants_a_num() const noexcept {return participants_a_num_;}
    std::size_t participants_b_num() const noexcept {return participants_b_.size();}

    real_type max_cutoff_length() const noexcept
    {
        return v0_range_;
    }

    void initialize(const system_type& sys, const topology_type& topol) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->update(sys, topol);
        is_following_participants_.resize(sys.size(), false);
        for(const auto& i : participants_b_)
        {
            is_following_participants_[i] = true;
        }
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
        return make_range(participants_.begin(), participants_.end());
    }
    range<typename std::vector<std::size_t>::const_iterator>
    possible_partners_of(const std::size_t /*participant_idx*/,
                         const std::size_t particle_idx) const noexcept
    {
        if(is_following_participants_[particle_idx])
        {
            return make_range(participants_a_.begin(), participants_a_.end());
        }
        else
        {
            return make_range(participants_b_.begin(), participants_b_.end());
        }
    }
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept
    {
        // i is leading_participant and j is a member of possible_partners_of
        // from implementation of SpatialPartition
        if(is_following_participants_[i])
        {
            if(!is_following_participants_[j])
            {
                return !exclusion_list_.is_excluded(i, j);
            }
        }
        else if(is_following_participants_[j])
        {
            return !exclusion_list_.is_excluded(i, j);
        }
        return false;
    }

    // -----------------------------------------------------------------------
    // used by Observer.
    static const char* name() noexcept {return "GlobalStoichiometricInteraction";}

  private:

    real_type                v0_;
    real_type                v0_range_;
    real_type                inv_range_;
    real_type                inv_range_6_;
    std::vector<std::size_t> participants_;
    std::vector<std::size_t> participants_a_;
    std::size_t              participants_a_num_;
    std::vector<std::size_t> participants_b_;
    std::vector<bool>        is_following_participants_;

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
