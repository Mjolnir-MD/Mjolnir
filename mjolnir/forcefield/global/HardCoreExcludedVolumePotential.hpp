#ifndef MJOLNIR_POTENTIAL_GLOBAL_HARD_CORE_EXCLUDED_VOLUME_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_HARD_CORE_EXCLUDED_VOLUME_POTENTIAL_HPP
#include <mjolnir/forcefield/global/ParameterList.hpp>
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>

namespace mjolnir
{

// hard core excluded volume potential.
// This class contains radii of hard core, thicknesses of soft layer represented
// by default excluded volume, calculates energy and derivative of the potential
// function.
template<typename realT>
class HardCoreExcludedVolumePotential
{
  public:
    using real_type      = realT;
    using parameter_type = std::tuple<real_type, real_type, real_type>; // {sigma, epsilon, radii}
    using self_type      = HardCoreExcludedVolumePotential<real_type>;

    static constexpr real_type default_cutoff() noexcept
    {
        return real_type(2.0);
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{real_type(0.0), real_type(0.0)};
    }

    static void set_cutoff_ratio(const real_type ratio)
    {
        if(self_type::cutoff_ratio < ratio)
        {
            self_type::cutoff_ratio   = ratio;
            self_type::coef_at_cutoff = std::pow(real_type(1) / ratio, 12);
        }
        return;
    }

    static real_type cutoff_ratio;
    static real_type coef_at_cutoff;

  public:

    explicit HardCoreExcludedVolumePotential(const parameter_type& params) noexcept
        : sigma_(std::get<0>(params)), epsilon_(std::get<1>(params)), radii_(std::get<2>(params))
    {}
    ~HardCoreExcludedVolumePotential() = default;
    HardCoreExcludedVolumePotential(const HardCoreExcludedVolumePotential&) = default;
    HardCoreExcludedVolumePotential(HardCoreExcludedVolumePotential&&)      = default;
    HardCoreExcludedVolumePotential& operator=(const HardCoreExcludedVolumePotential&) = default;
    HardCoreExcludedVolumePotential& operator=(HardCoreExcludedVolumePotential&&)      = default;

    real_type potential(const real_type r) const noexcept
    {
        if(this->sigma_ * self_type::cutoff_ratio + this->radii_ < r)
        {
            return 0;
        }
        else if(r < this->radii_)
        {
            return std::numeric_limits<real_type>::infinity();
        }

        const real_type gap = r - this->radii_;
        const real_type sigma1gap1 = this->sigma_ / gap;
        const real_type sigma3gap3 = sigma1gap1 * sigma1gap1 * sigma1gap1;
        const real_type sigma6gap6 = sigma3gap3 * sigma3gap3;

        return this->epsilon_ * (sigma6gap6 * sigma6gap6 - self_type::coef_at_cutoff);
    }
    real_type derivative(const real_type r) const noexcept
    {
        if(this->sigma_ * self_type::cutoff_ratio + this->radii_ < r)
        {
            return 0.0;
        }
        else if(r < this->radii_)
        {
            return std::numeric_limits<real_type>::infinity();
        }
        const real_type gapinv = 1 / (r - this->radii_);
        const real_type sigma1gap1 = this->sigma_ * gapinv;
        const real_type sigma3gap3 = sigma1gap1 * sigma1gap1 * sigma1gap1;
        const real_type sigma6gap6 = sigma3gap3 * sigma3gap3;

        return -12.0 * this->epsilon_ * sigma6gap6 * sigma6gap6 * gapinv;
    }

    template<typename T>
    void initialize(const System<T>&) noexcept {return;}

    template<typename T>
    void update(const System<T>&) noexcept {return;}

    real_type cutoff() const noexcept {return sigma_ * self_type::cutoff_ratio;}

    // ------------------------------------------------------------------------
    // used by Observer.

    static const char* name() noexcept {return "HardCoreExcludedVolume";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    real_type& epsilon()       noexcept {return this->epsilon_;}
    real_type  epsilon() const noexcept {return this->epsilon_;}
    real_type& sigma()         noexcept {return this->sigma_;}
    real_type  sigma()   const noexcept {return this->sigma_;}
    real_type& radii()         noexcept {return this->radii_;}
    real_type  radii()   const noexcept {return this->radii_;}

  private:

    real_type sigma_; // sigma_i + sigma_j
    real_type epsilon_;
    real_type radii_; // radius_i + radius_j
};
template<typename realT>
typename HardCoreExcludedVolumePotential<realT>::real_type
HardCoreExcludedVolumePotential<realT>::cutoff_ratio  = 0.0;
template<typename realT>
typename HardCoreExcludedVolumePotential<realT>::real_type
HardCoreExcludedVolumePotential<realT>::coef_at_cutoff = 0.0;


template<typename traitsT>
class HardCoreExcludedVolumeParameterList
    : public ParameterListBase<traitsT, HardCoreExcludedVolumePotential<typename traitsT::real_type>>
{
  public:
    using traits_type    = traitsT;
    using real_type      = typename traits_type::real_type;
    using potential_type = HardCoreExcludedVolumePotential<real_type>;
    using base_type      = ParameterListBase<traits_type, potential_type>;


    // `pair_parameter_type` is a parameter for a interacting pair.
    // Although it is the same type as `parameter_type` in this potential,
    // it can be different from normal parameter for each particle.
    using parameter_type      = std::pair<real_type, real_type>; // {sigma, hardcore_radius}
    using pair_parameter_type = typename potential_type::parameter_type;
    using container_type      = std::vector<parameter_type>;

    // topology stuff
    using system_type    = System<traits_type>;
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList <traits_type>;

    static constexpr real_type default_cutoff() noexcept
    {
        return real_type(2.0);
    }

    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{0.0, 0.0};
    }

  public:

    HardCoreExcludedVolumeParameterList(
        const real_type eps, const real_type cutoff_ratio,
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : epsilon_(eps), cutoff_ratio_(cutoff_ratio),
          exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {
        this->parameters_.reserve(parameters.size());
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
    ~HardCoreExcludedVolumeParameterList() = default;
    HardCoreExcludedVolumeParameterList(const HardCoreExcludedVolumeParameterList&) = default;
    HardCoreExcludedVolumeParameterList(HardCoreExcludedVolumeParameterList&&) = default;
    HardCoreExcludedVolumeParameterList& operator=(const HardCoreExcludedVolumeParameterList&) = default;
    HardCoreExcludedVolumeParameterList& operator=(HardCoreExcludedVolumeParameterList&&) = default;

    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept override
    {
        const auto sigma1           = this->parameters_[i].first;
        const auto hardcore_radius1 = this->parameters_[i].second;
        const auto sigma2           = this->parameters_[j].first;
        const auto hardcore_radius2 = this->parameters_[j].second;

        return std::make_tuple(sigma1 + sigma2, this->epsilon_,
                               hardcore_radius1 + hardcore_radius2);
    }

    void initialize(const system_type& sys, const topology_type& topol) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->update(sys, topol);
        return;
    }

    void update(const system_type& sys, const topology_type& topol) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        if(parameters_.empty())
        {
            // dummy value that does not screw up cell lists
            this->max_cutoff_length_ = 1.0;
        }
        else
        {
            const auto max_iter = std::max_element(
                parameters_.begin(), parameters_.end(),
                [this](const parameter_type& lhs, const parameter_type& rhs) noexcept {
                    const auto rc1 = cutoff_ratio_ * lhs.first + lhs.second;
                    const auto rc2 = cutoff_ratio_ * rhs.first + rhs.second;
                    return rc1 < rc2;
                });
            this->max_cutoff_length_ =
                cutoff_ratio_ * (max_iter->first * 2) + max_iter->second;
        }

        // update exclusion list based on sys.topology()
        exclusion_list_.make(sys, topol);
        return;
    }

    real_type cutoff_ratio()      const noexcept {return this->cutoff_ratio_;}
    real_type max_cutoff_length() const noexcept override {return this->max_cutoff_length_;}

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
        return make_range(participants_.begin() + participant_idx + 1,
                          participants_.end());
    }
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept override
    {
        return (i < j) && !exclusion_list_.is_excluded(i, j);
    }

    // for testing
    exclusion_list_type const& exclusion_list() const noexcept override
    {
        return exclusion_list_;
    }

    // ------------------------------------------------------------------------
    // used by Observer.

    static const char* name() noexcept {return "HardCoreExcludedVolume";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.
    real_type& epsilon()       noexcept {return this->epsilon_;}
    real_type  epsilon() const noexcept {return this->epsilon_;}
    std::vector<parameter_type>&       parameters()       noexcept {return parameters_;}
    std::vector<parameter_type> const& parameters() const noexcept {return parameters_;}

    base_type* clone() const override
    {
        return new HardCoreExcludedVolumeParameterList(*this);
    }

 private:

    real_type epsilon_;
    real_type cutoff_ratio_;
    real_type max_cutoff_length_;

    std::vector<parameter_type> parameters_;
    std::vector<std::size_t> participants_;

    exclusion_list_type exclusion_list_;
};
} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class HardCoreExcludedVolumeParameterList<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class HardCoreExcludedVolumeParameterList<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class HardCoreExcludedVolumeParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class HardCoreExcludedVolumeParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_HARD_CORE_EXCLUDED_VOLUME_POTENTIAL */
