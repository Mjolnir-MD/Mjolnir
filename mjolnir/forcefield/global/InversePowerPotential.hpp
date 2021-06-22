#ifndef MJOLNIR_POTENTIAL_GLOBAL_INVERSE_POWER_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_INVERSE_POWER_POTENTIAL_HPP
#include <mjolnir/forcefield/global/ParameterList.hpp>
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <algorithm>
#include <tuple>
#include <numeric>
#include <memory>
#include <cmath>

namespace mjolnir
{

// U = 1 / r^N
template<typename realT>
class InversePowerPotential
{
  public:

    using real_type    = realT;
    using integer_type = std::int32_t;

    struct parameter_type
    {
        real_type sigma;
    };

    static real_type default_cutoff(const integer_type n) noexcept
    {
        return std::pow(2.0, 12 / n);
    }

  public:

    InversePowerPotential(const real_type cutoff_ratio, const integer_type n,
                          const real_type epsilon)
        : n_(n), epsilon_(epsilon), cutoff_ratio_(cutoff_ratio),
          coef_at_cutoff_(std::pow(real_type(1) / cutoff_ratio, n))
    {}
    ~InversePowerPotential() = default;

    real_type potential(const real_type r, const parameter_type& params) const noexcept
    {
        if(params.sigma * cutoff_ratio_ < r) {return 0.0;}

        const real_type d_r = params.sigma / r;
        const real_type drn = std::pow(d_r, n_);
        return epsilon_ * (drn - coef_at_cutoff_);
    }
    real_type derivative(const real_type r, const parameter_type& params) const noexcept
    {
        if(params.sigma * cutoff_ratio_ < r) {return 0.0;}

        const real_type rinv = 1.0 / r;
        const real_type d_r = params.sigma * rinv;
        const real_type drn = std::pow(d_r, n_);
        return -n_ * epsilon_ * drn * rinv;
    }

    template<typename T>
    void initialize(const System<T>&) noexcept {return;}

    template<typename T>
    void update(const System<T>&) noexcept {return;}

    // ------------------------------------------------------------------------
    // used by Observer.

    static const char* name() noexcept {return "InversePower";}

    // It takes per-particle parameters and return the maximum cutoff length.
    // CombinationRule normally uses this.
    // Note that, pair-parameter and per-particle parameter can differ from
    // each other. Lorentz-Bertherot uses the same parameter_type because it is
    // for L-J and L-J-like potentials that has {sigma, epsilon} for each
    // particle and also for each pair of particles.
    template<typename InputIterator>
    real_type max_cutoff(const InputIterator first, const InputIterator last) const noexcept
    {
        static_assert(std::is_same<
                typename std::iterator_traits<InputIterator>::value_type,
                parameter_type>::value, "");

        if(first == last) {return 1;}

        real_type max_sigma = 0;
        for(auto iter = first; iter != last; ++iter)
        {
            const auto& parameter = *iter;
            max_sigma = std::max(max_sigma, parameter.sigma);
        }
        return max_sigma * cutoff_ratio_;
    }
    // It returns absolute cutoff length using pair-parameter.
    // `CombinationTable` uses this.
    real_type absolute_cutoff(const parameter_type& params) const noexcept
    {
        return params.sigma * cutoff_ratio_;
    }

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    real_type& epsilon()       noexcept {return this->epsilon_;}
    real_type  epsilon() const noexcept {return this->epsilon_;}

    integer_type& n()       noexcept {return this->n_;}
    integer_type  n() const noexcept {return this->n_;}

    real_type  cutoff_ratio() const noexcept {return this->cutoff_ratio_;}

  private:

    integer_type n_;
    real_type epsilon_;
    real_type cutoff_ratio_;
    real_type coef_at_cutoff_;
};

template<typename traitsT>
class InversePowerParameterList
    : public ParameterListBase<traitsT, InversePowerPotential<typename traitsT::real_type>>
{
  public:
    using traits_type    = traitsT;
    using real_type      = typename traits_type::real_type;
    using potential_type = InversePowerPotential<real_type>;
    using base_type      = ParameterListBase<traits_type, potential_type>;
    using self_type      = InversePowerParameterList<traits_type>;

    // `pair_parameter_type` is a parameter for an interacting pair.
    // Although it is the same type as `parameter_type` in this potential,
    // it can be different from normal parameter for each particle.
    // This enables NeighborList to cache a value to calculate forces between
    // the particles, e.g. by having qi * qj for pair of particles i and j.

    using parameter_type       = typename potential_type::parameter_type;
    using pair_parameter_type  = typename potential_type::parameter_type;
    using container_type       = std::vector<parameter_type>;

    static constexpr parameter_type default_parameter()
    {
        return parameter_type{real_type(0)};
    }

    // topology stuff
    using system_type          = System<traits_type>;
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList<traits_type>;

  public:

    InversePowerParameterList(
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>&         exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
      : exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {
        this->parameters_  .reserve(parameters.size());
        this->participants_.reserve(parameters.size());
        for(const auto& idxp : parameters)
        {
            const auto idx = idxp.first;
            this->participants_.push_back(idx);
            if(idx >= this->parameters_.size())
            {
                this->parameters_.resize(idx+1, self_type::default_parameter());
            }
            this->parameters_.at(idx) = idxp.second;
        }
    }
    ~InversePowerParameterList() = default;
    InversePowerParameterList(const InversePowerParameterList&) = default;
    InversePowerParameterList(InversePowerParameterList&&)      = default;
    InversePowerParameterList& operator=(const InversePowerParameterList&) = default;
    InversePowerParameterList& operator=(InversePowerParameterList&&)      = default;

    // this value will be stored in NeighborList
    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept override
    {
        return pair_parameter_type{parameters_[i].sigma + parameters_[j].sigma};
    }

    void initialize(const system_type& sys, const topology_type& topol,
                    const potential_type& pot) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->update(sys, topol, pot);
        return;
    }

    // nothing to be done if system parameter (e.g. temperature) do not changes
    void update(const system_type& sys, const topology_type& topol,
                const potential_type& pot) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->max_cutoff_length_ = pot.max_cutoff(parameters_.begin(), parameters_.end());
        exclusion_list_.make(sys, topol);
        return;
    }

    real_type max_cutoff_length() const noexcept override
    {
        return max_cutoff_length_;
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

    std::vector<std::size_t> const& participants() const noexcept override {return participants_;}

    range<typename std::vector<std::size_t>::const_iterator>
    leading_participants() const noexcept override
    {
        return make_range(participants_.begin(), std::prev(participants_.end()));
    }
    range<typename std::vector<std::size_t>::const_iterator>
    possible_partners_of(const std::size_t participant_idx,
                         const std::size_t /*particle_idx*/) const noexcept override
    {
        return make_range(participants_.begin() + participant_idx + 1,
                          participants_.end());
    }
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept override
    {
        return (i < j) && !exclusion_list_.is_excluded(i, j);
    }

    // ------------------------------------------------------------------------
    // for testing
    exclusion_list_type const& exclusion_list() const noexcept override
    {
        return exclusion_list_;
    }

    // ------------------------------------------------------------------------
    // used by Observer.

    static const char* name() noexcept {return "InversePower";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    std::vector<parameter_type>&       parameters()       noexcept {return this->parameters_;}
    std::vector<parameter_type> const& parameters() const noexcept {return this->parameters_;}

    base_type* clone() const override
    {
        return new InversePowerParameterList(*this);
    }

  private:

    real_type max_cutoff_length_;
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
extern template class InversePowerParameterList<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class InversePowerParameterList<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class InversePowerParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class InversePowerParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_EXCLUDED_VOLUME_POTENTIAL */
