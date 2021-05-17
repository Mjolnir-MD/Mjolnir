#ifndef MJOLNIR_FORCEFIELD_GLOBAL_TABULATED_LENNARD_JONES_ATTRACTIVE_POTENTIAL_HPP
#define MJOLNIR_FORCEFIELD_GLOBAL_TABULATED_LENNARD_JONES_ATTRACTIVE_POTENTIAL_HPP
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/logger.hpp>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace mjolnir
{

template<typename traitsT>
class TabulatedLennardJonesAttractivePotential
{
  public:
    using traits_type          = traitsT;
    using real_type            = typename traits_type::real_type;
    using system_type          = System<traits_type>;
    using parameter_type       = std::string;
    using pair_parameter_type  = std::pair<real_type, real_type>; // {sigma, epsilon};
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
        return real_type(2.5);
    }
    static parameter_type default_parameter() noexcept
    {
        return std::string("");
    }

  public:

    TabulatedLennardJonesAttractivePotential(
        const real_type cutoff_ratio, table_type&& table,
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>&         exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
    : cutoff_ratio_(cutoff_ratio),
      coef_at_cutoff_(std::pow(1 / cutoff_ratio, 12) - std::pow(1 / cutoff_ratio, 6)),
      table_(std::move(table)),
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
    ~TabulatedLennardJonesAttractivePotential() = default;

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
        constexpr real_type sixth_root_of_two(1.12246204831);

        const real_type sigma   = p.first;
        const real_type epsilon = p.second;
        if(r < sigma * sixth_root_of_two)
        {
            return -epsilon;
        }
        else if(sigma * this->cutoff_ratio_ < r)
        {
            return 0;
        }

        const real_type sr1 = sigma / r;
        const real_type sr3 = sr1 * sr1 * sr1;
        const real_type sr6 = sr3 * sr3;
        return 4 * epsilon * (sr6 * (sr6 - 1) - coef_at_cutoff_);
    }
    real_type derivative(const real_type r, const pair_parameter_type& p) const noexcept
    {
        constexpr real_type sixth_root_of_two(1.12246204831);

        const real_type sigma   = p.first;
        const real_type epsilon = p.second;
        if(r < sigma * sixth_root_of_two ||
               sigma * this->cutoff_ratio_ < r)
        {
            return 0;
        }

        const real_type rinv = 1 / r;
        const real_type sr1 = sigma * rinv;
        const real_type sr3 = sr1 * sr1 * sr1;
        const real_type sr6 = sr3 * sr3;
        return 24 * epsilon * (sr6 - 2 * sr6 * sr6) * rinv;
    }

    real_type cutoff_ratio()   const noexcept {return this->cutoff_ratio_;}
    real_type coef_at_cutoff() const noexcept {return this->coef_at_cutoff_;}

    real_type max_cutoff_length() const noexcept
    {
        using value_type = typename table_type::value_type;
        const real_type max_sigma = std::max_element(
            this->table_.cbegin(), this->table_.cend(),
            [](const value_type& lhs, const value_type& rhs) noexcept {
                return lhs.second.first < rhs.second.first;
            })->second.first;

        return max_sigma * this->cutoff_ratio_;
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
    static const char* name() noexcept {return "TabulatedLennardJonesAttractive";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    // access to the parameters...
    std::vector<parameter_type>&       parameters()       noexcept {return parameters_;}
    std::vector<parameter_type> const& parameters() const noexcept {return parameters_;}

    table_type&       table()       noexcept {return table_;}
    table_type const& table() const noexcept {return table_;}

  private:

    real_type cutoff_ratio_;
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
extern template class TabulatedLennardJonesAttractivePotential<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class TabulatedLennardJonesAttractivePotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class TabulatedLennardJonesAttractivePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class TabulatedLennardJonesAttractivePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
