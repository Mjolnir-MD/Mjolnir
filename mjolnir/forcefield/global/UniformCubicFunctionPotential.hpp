#ifndef MJOLNIR_POTENTIAL_GLOBAL_CUBIC_FUNCTION_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_CUBIC_FUNCTION_POTENTIAL_HPP
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/util/empty.hpp>
#include <mjolnir/util/logger.hpp>
#include <vector>

namespace mjolnir
{

// This potential is same to GlobalStoichiometricInteractionPotential
// and exist to use that in non-stoichiometric interaction.
template<typename traitsT>
class UniformCubicFunctionPotential
{
  public:
    using traits_type         = traitsT;
    using real_type           = typename traits_type::real_type;
    using system_type         = System<traits_type>;
    using parameter_type      = empty_t; // no particle_specific parameter
    using container_type      = empty_t; // no parameter, so no container there.
    using pair_parameter_type = empty_t; // no particle-pair-specific parameter

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
    UniformCubicFunctionPotential(
        const real_type epsilon, const real_type v0, const real_type range,
        const std::vector<std::pair<std::size_t, parameter_type>>&   parameters,
        const std::map<connection_kind_type, std::size_t>&           exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : epsilon_(epsilon), v0_(v0), range_(range), v0_range_(v0 + range),
          inv_range_(1.0 / range), inv_range_6_(6.0 / range),
          exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {
        MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
        MJOLNIR_LOG_FUNCTION_DEBUG();

        this->participants_.reserve(parameters.size());
        for(const auto& idxp : parameters)
        {
            this->participants_.push_back(idxp.first);
        }
    }
    ~UniformCubicFunctionPotential() = default;

    pair_parameter_type prepare_params(std::size_t, std::size_t) const noexcept
    {
        return pair_parameter_type{}; // no pre-calculated parameter
    }

    // forwarding functions for clarity...
    real_type potential(const std::size_t i, const std::size_t j,
                        const real_type r) const noexcept
    {
        return this->potential(r, this->prepare_params(i, j));
    }
    real_type derivative(const std::size_t i, const std::size_t j,
                         const real_type r) const noexcept
    {
        return this->derivative(r, this->prepare_params(i, j));
    }

    real_type potential(const real_type r, const pair_parameter_type&) const noexcept
    {
        if(r < v0_){return -1.0;}
        else if(v0_range_ < r){return 0.0;}

        const real_type r_v0           = r - v0_;
        const real_type r_v0_inv_range = r_v0 * inv_range_;
        const real_type r_v0_inv_range_2  = r_v0_inv_range * r_v0_inv_range;
        const real_type r_v0_inv_range_3  = r_v0_inv_range_2 * r_v0_inv_range;
        return -epsilon_ * (1.0 - 3.0 * r_v0_inv_range_2 + 2.0 * r_v0_inv_range_3);
    }

    real_type derivative(const real_type r, const pair_parameter_type&) const noexcept
    {
        if(r < v0_ || v0_range_ < r){return 0.0;}

        const real_type r_v0          = r - v0_;
        const real_type r_v0_inv_range   = r_v0 * inv_range_;
        const real_type r_v0_inv_range_2 = r_v0_inv_range * r_v0_inv_range;
        return -epsilon_ * inv_range_6_ * (r_v0_inv_range_2 - r_v0_inv_range);
    }

    real_type max_cutoff_length() const noexcept
    {
        return v0_range_;
    }

    void initialize(const system_type& sys, const topology_type& topol) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // if no participants are given, consider all the particles are related.
        if(this->participants_.empty())
        {
            MJOLNIR_LOG_WARN("UniformCubicFunctionPotential does not have any participants.");
            MJOLNIR_LOG_WARN("All the particles are considered to be participants.");
            MJOLNIR_LOG_WARN("To disable this potential, "
                             "simply remove the part from the input file.");

            this->participants_.resize(sys.size());
            std::iota(this->participants_.begin(), this->participants_.end(), 0u);
        }

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
        return make_range(participants_.begin(), participants_.end());
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
        // if not excluded, the pair has interaction.
        return (i < j) && !exclusion_list_.is_excluded(i, j);
    }
    exclusion_list_type const& exclusion_list() const noexcept
    {
        return exclusion_list_;
    }

    // -----------------------------------------------------------------------
    // used by Observer.
    static const char* name() noexcept {return "UniformCubicFunction";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    // access to the parameters...
    real_type& epsilon()                   noexcept {return epsilon_;}
    real_type  epsilon()             const noexcept {return epsilon_;}
    real_type& v0()                        noexcept {return v0_;}
    real_type  v0()                  const noexcept {return v0_;}
    real_type& interaction_range()         noexcept {return range_;}
    real_type  interaction_range()   const noexcept {return range_;}

  private:

    real_type                epsilon_;
    real_type                v0_;
    real_type                range_;
    real_type                v0_range_;
    real_type                inv_range_;
    real_type                inv_range_6_;
    std::vector<std::size_t> participants_;

    exclusion_list_type exclusion_list_;
};
} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class UniformCubicFunctionPotential<SimulatorTraits<double, UnlimitedBoundary       >>;
extern template class UniformCubicFunctionPotential<SimulatorTraits<float,  UnlimitedBoundary       >>;
extern template class UniformCubicFunctionPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class UniformCubicFunctionPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_POTENTIAL_GLOBAL_CUBIC_FUNCTION_POTENTIAL_HPP */
