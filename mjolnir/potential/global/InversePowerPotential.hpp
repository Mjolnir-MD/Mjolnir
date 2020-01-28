#ifndef MJOLNIR_POTENTIAL_GLOBAL_INVERSE_POWER_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_INVERSE_POWER_POTENTIAL_HPP
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <algorithm>
#include <numeric>
#include <memory>
#include <cmath>

namespace mjolnir
{

// inverse power potential.
template<typename traitsT>
class InversePowerPotential
{
  public:

    using traits_type           = traitsT;
    using real_type             = typename traits_type::real_type;
    using integer_type          = std::int32_t;
    using parameter_type        = real_type;
    using container_type        = std::vector<parameter_type>;

    using pair_parameter_type   = parameter_type;

    // ------------------------------------------------------------------------
    // topology stuff
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList <traits_type>;

    // std::pow is not marked constexpr
    static real_type default_cutoff(const integer_type n) noexcept
    {
        return std::pow(2.0, 12.0 / n);
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{0.0};
    }

  public:

    InversePowerPotential(const real_type eps, const integer_type n,
        const real_type cutoff_ratio,
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>&         exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        :epsilon_(eps), n_(n), cutoff_ratio_(cutoff_ratio),
         coef_at_cutoff_(std::pow(1.0 / cutoff_ratio, n)),
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
    ~InversePowerPotential() = default;
    InversePowerPotential(const InversePowerPotential&) = default;
    InversePowerPotential(InversePowerPotential&&)      = default;
    InversePowerPotential& operator=(const InversePowerPotential&) = default;
    InversePowerPotential& operator=(InversePowerPotential&&)      = default;

    // this value will be stored in NeighborList
    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept
    {
        return this->parameters_[i] + this->parameters_[j];
    }

    // forwarding function for clarity
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

    real_type potential(const real_type r, const pair_parameter_type& d) const noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
        MJOLNIR_LOG_FUNCTION_DEBUG();

        if(d * this->cutoff_ratio_ < r){return 0.0;}

        const real_type d_r = d / r;
        const real_type drn = std::pow(d_r, n_);
        return this->epsilon_ * (drn - this->coef_at_cutoff_);
    }
    real_type derivative(const real_type r, const pair_parameter_type& d) const noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
        MJOLNIR_LOG_FUNCTION_DEBUG();

        if(d * this->cutoff_ratio_ < r){return 0.0;}

        const real_type rinv = 1.0 / r;
        const real_type d_r = d * rinv;
        const real_type drn = std::pow(d_r, n_);
        return -n_ * this->epsilon_ * drn * rinv;
    }

    void initialize(const System<traits_type>& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->update(sys);
        return;
    }

    // nothing to be done if system parameter (e.g. temperature) do not changes
    void update(const System<traits_type>& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // update exclusion list based on sys.topology()
        exclusion_list_.make(sys, sys.topology());
        return;
    }

    real_type cutoff_ratio()   const noexcept {return this->cutoff_ratio_;}
    real_type coef_at_cutoff() const noexcept {return this->coef_at_cutoff_;}

    real_type max_cutoff_length() const
    {
        const real_type max_sigma =
            * (std::max_element(parameters_.cbegin(), parameters_.cend()));
        return 2 * max_sigma * this->cutoff_ratio_;
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

    // ------------------------------------------------------------------------
    // for testing
    exclusion_list_type const& exclusion_list() const noexcept
    {
        return exclusion_list_;
    }

    // ------------------------------------------------------------------------
    // used by Observer.

    static const char* name() noexcept {return "InversePower";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    real_type& epsilon()       noexcept {return this->epsilon_;}
    real_type  epsilon() const noexcept {return this->epsilon_;}
    integer_type& n()             noexcept {return this->n_;}
    integer_type  n()       const noexcept {return this->n_;}
    std::vector<real_type>&       parameters()       noexcept {return this->parameters_;}
    std::vector<real_type> const& parameters() const noexcept {return this->parameters_;}


  private:

    real_type epsilon_;
    integer_type n_;
    real_type cutoff_ratio_;
    real_type coef_at_cutoff_; // correction of energy value at cutoff length
    std::vector<parameter_type> parameters_;
    std::vector<std::size_t>    participants_;

    exclusion_list_type exclusion_list_;
};
} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class InversePowerPotential<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class InversePowerPotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class InversePowerPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class InversePowerPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_EXCLUDED_VOLUME_POTENTIAL */
