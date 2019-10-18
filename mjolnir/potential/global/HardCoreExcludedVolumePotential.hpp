#ifndef MJOLNIR_POTENTIAL_GLOBAL_HARD_CORE_EXCLUDED_VOLUME_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_HARD_CORE_EXCLUDED_VOLUME_POTENTIAL_HPP
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
template<typename traitsT>
class HardCoreExcludedVolumePotential
{
  public:
    using traits_type    = traitsT;
    using real_type      = typename traits_type::real_type;
    using parameter_type = std::pair<real_type, real_type>; // {sigma, hardcore_radius}
    using container_type = std::vector<parameter_type>;
    // `pair_parameter_type` is a parameter for a interacting pair.
    // Although it is the same type as `parameter_type` in this potential,
    // it can be different from normal parameter for each particle.
    using pair_parameter_type = parameter_type;

    // topology stuff
    using topology_type = Topology;
    using molecule_id_type = typename topology_type::molecule_id_type;
    using group_id_type = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList;

    static constexpr real_type default_cutoff() noexcept
    {
        return real_type(2.0);
    }

    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{0.0, 0.0};
    }

  public:

    HardCoreExcludedVolumePotential(const real_type eps, const real_type cutoff_ratio,
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : epsilon_(eps), cutoff_ratio_(cutoff_ratio),
          coef_at_cutoff_(std::pow(1.0 / cutoff_ratio, 12)),
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
    ~HardCoreExcludedVolumePotential() = default;
    HardCoreExcludedVolumePotential(const HardCoreExcludedVolumePotential&) = default;
    HardCoreExcludedVolumePotential(HardCoreExcludedVolumePotential&&) = default;
    HardCoreExcludedVolumePotential& operator=(const HardCoreExcludedVolumePotential&) = default;
    HardCoreExcludedVolumePotential& operator=(HardCoreExcludedVolumePotential&&) = default;

    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept
    {
        const auto sigma1           = this->parameters_[i].first;
        const auto hardcore_radius1 = this->parameters_[i].second;
        const auto sigma2           = this->parameters_[j].first;
        const auto hardcore_radius2 = this->parameters_[j].second;

        return std::make_pair(sigma1 + sigma2, hardcore_radius1 + hardcore_radius2);
    }

    //forwarding fumctions for clarity...
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

    real_type potential(const real_type r, const pair_parameter_type& pram) const noexcept
    {
        const real_type sigma_sum = pram.first;
        const real_type hradi_sum = pram.second;
        if(sigma_sum * this->cutoff_ratio_ + hradi_sum < r){return 0;}
        else if(r < hradi_sum){return std::numeric_limits<real_type>::infinity();}

        const real_type gap = r - hradi_sum;
        const real_type sigma1gap1 = sigma_sum / gap;
        const real_type sigma3gap3 = sigma1gap1 * sigma1gap1 * sigma1gap1;
        const real_type sigma6gap6 = sigma3gap3 * sigma3gap3;
        const real_type sigma12gap12 = sigma6gap6 * sigma6gap6;

        return this->epsilon_ * (sigma12gap12 - this-> coef_at_cutoff_);
    }
    real_type derivative(const real_type r, const pair_parameter_type& pram) const noexcept
    {
        const real_type sigma_sum = pram.first;
        const real_type hradi_sum = pram.second;
        if(sigma_sum * this->cutoff_ratio_ + hradi_sum < r){return 0.0;}
        else if(r < hradi_sum){return std::numeric_limits<real_type>::infinity();}
        const real_type gapinv = 1 / (r - hradi_sum);
        const real_type sigma1gap1 = sigma_sum * gapinv;
        const real_type sigma3gap3 = sigma1gap1 * sigma1gap1 * sigma1gap1;
        const real_type sigma6gap6 = sigma3gap3 * sigma3gap3;
        const real_type sigma12gap12 = sigma6gap6 * sigma6gap6;
        return -12.0 * this->epsilon_ * sigma12gap12 * gapinv;
    }

    void initialize(const System<traits_type>& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->update(sys);
        return;
    }

    void update(const System<traits_type>& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // update exclusion list based on sys.topology()
        exclusion_list_.make(sys);
    }

    real_type cutoff_ratio() const noexcept {return this->cutoff_ratio_;}
    real_type coef_at_cutoff() const noexcept {return this->coef_at_cutoff_;}

    real_type max_cutoff_length() const noexcept
    {
        const parameter_type& max_param = *std::max_element(
            this->parameters_.cbegin(), this->parameters_.cend(),
            [this](const parameter_type& lhs, const parameter_type& rhs) noexcept {
                return this->cutoff_length(lhs) < this->cutoff_length(rhs);
            });
        return max_param.first * this->cutoff_ratio_ + max_param.second;
    }

    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept
    {
        // if not excluded, the pair has interaction.
        return !exclusion_list_.is_excluded(i, j);
    }
    // for testing
    exclusion_list_type const& exclusion_list() const noexcept
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

    std::vector<std::size_t> const& participants() const noexcept {return participants_;}

 private:
    real_type epsilon_;
    real_type hard_raius;
    real_type cutoff_ratio_;
    real_type coef_at_cutoff_;
    std::vector<parameter_type> parameters_;
    std::vector<std::size_t> participants_;

    exclusion_list_type exclusion_list_;

    real_type cutoff_length(const parameter_type& pram) const noexcept
    {
        return this->cutoff_ratio_ * pram.first + pram.second;
    }
};
} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class HardCoreExcludedVolumePotential<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class HardCoreExcludedVolumePotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class HardCoreExcludedVolumePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class HardCoreExcludedVolumePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_HARD_CORE_EXCLUDED_VOLUME_POTENTIAL */
