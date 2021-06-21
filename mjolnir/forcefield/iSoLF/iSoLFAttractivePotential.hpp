#ifndef MJOLNIR_POTENTIAL_GLOBAL_ISOLF_ATTRACTIVE_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_ISOLF_ATTRACTIVE_POTENTIAL_HPP
#include <mjolnir/forcefield/global/ParameterList.hpp>
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/math.hpp>
#include <vector>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <cmath>

namespace mjolnir
{

// attractive part of the iSoLF potential for coarse-grained lipids developed by
// - Diego Ugarte La Torre and Shoji Takada (2020) J. Chem. Phys 153, 205101
//   https://doi.org/10.1063/5.0026342
//
template<typename realT>
class iSoLFAttractivePotential
{
  public:
    using real_type = realT;

    struct parameter_type
    {
        real_type sigma;
        real_type epsilon;
        real_type omega;
        real_type romega;
    };

    static constexpr real_type default_cutoff() noexcept
    {
        return std::numeric_limits<real_type>::infinity();
    }

  public:

    iSoLFAttractivePotential() noexcept {}
    ~iSoLFAttractivePotential() = default;

    real_type potential(const real_type r, const parameter_type& params) const noexcept
    {
        constexpr real_type rc = 1.12246204831;
        constexpr real_type pi = math::constants<real_type>::pi();

        const real_type r_sigma_rc = r - params.sigma * rc; // r - sqrt[6]{2} sigma

        if     (r_sigma_rc <= 0)          {return -params.epsilon;}
        else if(params.omega < r_sigma_rc){return 0;}

        const real_type cosine = std::cos(pi * params.romega * r_sigma_rc * real_type(0.5));

        return -params.epsilon * cosine * cosine;
    }
    real_type derivative(const real_type r, const parameter_type& params) const noexcept
    {
        constexpr real_type rc = 1.12246204831;
        constexpr real_type pi = math::constants<real_type>::pi();

        const real_type r_sigma_rc = r - params.sigma * rc; // r - sqrt[6]{2} sigma

        if (r_sigma_rc <= 0 || params.omega < r_sigma_rc) {return 0;}

        const real_type sine = std::sin(pi * params.romega * r_sigma_rc);

        return params.epsilon * pi * params.romega * sine * real_type(0.5);
    }

    template<typename T>
    void initialize(const System<T>&) noexcept {return;}

    template<typename T>
    void update(const System<T>&) noexcept {return;}

    // ------------------------------------------------------------------------

    // It takes per-particle parameters and return the maximum cutoff length.
    // CombinationRule normally uses this.
    // Note that, pair-parameter contains romega to avoid extraneous division,
    // but the per-particle parameters does not have romega because it can be
    // derived from omega.
    template<typename InputIterator>
    real_type max_cutoff(const InputIterator first, const InputIterator last) const noexcept
    {
        if(first == last) {return 1;}
        constexpr real_type rc = 1.12246204831; // sqrt[6]{2}

        real_type max_sigma = 0;
        real_type max_omega = 0;
        for(auto iter = first; iter != last; ++iter)
        {
            const auto& parameter = *iter;
            max_sigma = std::max(max_sigma, parameter.sigma);
            max_omega = std::max(max_omega, parameter.omega);
        }
        return max_sigma * rc + max_omega;
    }
    // It returns absolute cutoff length using pair-parameter.
    // `CombinationTable` uses this.
    real_type absolute_cutoff(const parameter_type& params) const noexcept
    {
        constexpr real_type rc = 1.12246204831; // sqrt[6]{2}
        return params.sigma * rc + params.omega;
    }

    static const char* name() noexcept {return "iSoLFAttractive";}

    real_type cutoff_ratio()   const noexcept {return default_cutoff();}
    real_type coef_at_cutoff() const noexcept {return 0.0;}
};

template<typename traitsT>
class iSoLFAttractiveParameterList final
    : public ParameterListBase<traitsT, iSoLFAttractivePotential<typename traitsT::real_type>>
{
  public:
    using traits_type          = traitsT;
    using real_type            = typename traits_type::real_type;
    using potential_type       = iSoLFAttractivePotential<real_type>;
    using base_type            = ParameterListBase<traits_type, potential_type>;

    struct parameter_type
    {
        real_type sigma;
        real_type epsilon;
        real_type omega;
    };
    using pair_parameter_type  = typename potential_type::parameter_type;
    using container_type       = std::vector<parameter_type>;

    // topology stuff
    using system_type          = System<traits_type>;
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList <traits_type>;

    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{real_type(0), real_type(0), real_type(0)};
    }

  public:

    iSoLFAttractiveParameterList(
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
                this->parameters_.resize(idx+1, default_parameter());
            }
            this->parameters_.at(idx) = idxp.second;
        }
    }
    ~iSoLFAttractiveParameterList() = default;

    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept override
    {
        const auto sgm1 = parameters_[i].sigma;
        const auto eps1 = parameters_[i].epsilon;
        const auto omg1 = parameters_[i].omega;

        const auto sgm2 = parameters_[j].sigma;
        const auto eps2 = parameters_[j].epsilon;
        const auto omg2 = parameters_[j].omega;

        return pair_parameter_type{(sgm1 + sgm2) / 2,
                              ((eps1 == eps2) ? eps1 : std::sqrt(eps1 * eps2)),
                               (omg1 + omg2) / 2, real_type(2) / (omg1 + omg2)};
    }
    real_type max_cutoff_length() const noexcept override
    {
        return this->max_cutoff_length_;
    }

    void initialize(const system_type& sys, const topology_type& topol,
                    const potential_type& pot) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->update(sys, topol, pot);
        return;
    }

    void update(const system_type& sys, const topology_type& topol,
                const potential_type& pot) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->max_cutoff_length_ = pot.max_cutoff(parameters_.begin(), parameters_.end());
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

    std::vector<std::size_t> const& participants() const noexcept override
    {
        return participants_;
    }

    range<typename std::vector<std::size_t>::const_iterator>
    leading_participants() const noexcept override
    {
        return make_range(participants_.begin(), std::prev(participants_.end()));
    }
    range<typename std::vector<std::size_t>::const_iterator>
    possible_partners_of(const std::size_t participant_idx,
                         const std::size_t /*particle_idx*/) const noexcept override
    {
        return make_range(participants_.begin() + participant_idx + 1, participants_.end());
    }
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept override
    {
        // if not excluded, the pair has interaction.
        return (i < j) && !exclusion_list_.is_excluded(i, j);
    }
    // for testing
    exclusion_list_type const& exclusion_list() const noexcept override
    {
        return exclusion_list_;
    }

    // ------------------------------------------------------------------------
    // used by Observer.
    static const char* name() noexcept {return "iSoLFAttractive";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    // access to the parameters...
    std::vector<parameter_type>&       parameters()       noexcept {return parameters_;}
    std::vector<parameter_type> const& parameters() const noexcept {return parameters_;}

    base_type* clone() const override
    {
        return new iSoLFAttractiveParameterList(*this);
    }

  private:

    real_type max_cutoff_length_;
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
extern template class iSoLFAttractiveParameterList<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class iSoLFAttractiveParameterList<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class iSoLFAttractiveParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class iSoLFAttractiveParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_ISOLF_ATTRACTIVE_POTENTIAL */
