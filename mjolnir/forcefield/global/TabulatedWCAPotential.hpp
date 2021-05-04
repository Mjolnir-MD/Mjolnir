#ifndef MJOLNIR_POTENTIAL_GLOBAL_TABULATED_WCA_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_TABULATED_WCA_POTENTIAL_HPP
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/logger.hpp>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace mjolnir
{

// Well-known WCA interaction with Lorentz-Berthelot combining rules.
// This class contains sigmas and epsilons of the particles and calculates
// energy and derivative of the potential function.
template<typename traitsT>
class TabulatedWCAPotential
{
  public:
    using traits_type          = traitsT;
    using real_type            = typename traits_type::real_type;
    using system_type          = System<traits_type>;
    using parameter_type       = std::string;
    using pair_parameter_type  = std::pair<real_type, real_type>; // {sigma, epsilon}
    using container_type       = std::vector<parameter_type>;
    using table_type           = std::unordered_map<std::string, pair_parameter_type>;

    // topology stuff
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList <traits_type>;

    static constexpr real_type default_cutoff() noexcept
    {
        return real_type(1.12246204831); // pow(2.0, 1.0 / 6.0)
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return std::string("");
    }

  public:

    TabulatedWCAPotential(const real_type /*cutoff_ratio*/, table_type&& table,
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>&         exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
      : table_(std::move(table)),
        exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {
        this->parameters_  .reserve(parameters.size());
        this->participants_.reserve(parameters.size());
        for(const auto& idxp : parameters)
        {
            const auto idx = idxp.first;
            this->participants_.push_back(idx);
            if(idx >= this->parameters_.size())
            {
                this->parameters_.resize(idx+1, default_parameter());
            }
            this->parameters_.at(idx) = idxp.second;
        }
    }
    ~TabulatedWCAPotential() = default;

    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept
    {
        const auto key = parameters_[i] + ":" + parameters_[j];
        if(table_.count(key) == 0)
        {
            MJOLNIR_GET_DEFAULT_LOGGER();
            MJOLNIR_LOG_ERROR("parameter \"", key, "\" is not in the table");
        }
        return table_.at(key);
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

    real_type potential(const real_type r, const pair_parameter_type& p) const noexcept
    {
        const real_type sigma = p.first;
        if(sigma * default_cutoff() < r){return 0;}

        const real_type epsilon = p.second;

        const real_type r1s1   = sigma / r;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return 4 * epsilon * (r12s12 - r6s6 + real_type(0.25));
    }
    real_type derivative(const real_type r, const pair_parameter_type& p) const noexcept
    {
        const real_type sigma = p.first;
        if(sigma * default_cutoff() < r){return 0;}

        const real_type epsilon = p.second;

        const real_type rinv   = 1 / r;
        const real_type r1s1   = sigma * rinv;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return 24 * epsilon * (r6s6 - 2 * r12s12) * rinv;
    }

    real_type cutoff_ratio()   const noexcept {return default_cutoff();}
    real_type coef_at_cutoff() const noexcept {return 0.0;}

    real_type max_cutoff_length() const noexcept
    {
        using value_type = typename table_type::value_type;
        const real_type max_sigma = std::max_element(
            this->table_.cbegin(), this->table_.cend(),
            [](const value_type& lhs, const value_type& rhs) noexcept {
                return lhs.second.first < rhs.second.first;
            })->second.first;
        return max_sigma * this->default_cutoff();
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
        return;
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
        return make_range(participants_.begin() + participant_idx + 1, participants_.end());
    }
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept
    {
        // if not excluded, the pair has interaction.
        return (i < j) && !exclusion_list_.is_excluded(i, j);
    }
    // for testing
    exclusion_list_type const& exclusion_list() const noexcept
    {
        return exclusion_list_;
    }

    // ------------------------------------------------------------------------
    // used by Observer.
    static const char* name() noexcept {return "TabulatedWCA";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    // access to the parameters...
    std::vector<parameter_type>&       parameters()       noexcept {return parameters_;}
    std::vector<parameter_type> const& parameters() const noexcept {return parameters_;}

  private:

    real_type coef_at_cutoff_;
    table_type table_;
    container_type parameters_;
    std::vector<std::size_t> participants_;

    exclusion_list_type  exclusion_list_;
};

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class TabulatedWCAPotential<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class TabulatedWCAPotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class TabulatedWCAPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class TabulatedWCAPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
