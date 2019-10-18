#ifndef MJOLNIR_POTENTIAL_GLOBAL_EXCLUDED_VOLUME_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_EXCLUDED_VOLUME_POTENTIAL_HPP
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/math.hpp>
#include <algorithm>
#include <numeric>
#include <memory>
#include <cmath>

namespace mjolnir
{

// excluded volume potential.
// This class contains radii of the particles and calculates energy and
// derivative of the potential function.
// This class is an implementation of the excluded volume term used in
// Clementi's off-lattice Go-like model (Clement et al., 2000) and AICG2+ model
// (Li et al., 2014)
//
// Note: When ExcludedVolume is used with GlobalPairInteraction, `calc_force`
//       and `calc_energy` implemented here will not be used because we can
//       optimize the runtime efficiency by specializing GlobalPairInteraction.
//       See mjolnir/interaction/GlobalExcludedVolumeInteraction.hpp for detail.
template<typename realT>
class ExcludedVolumePotential
{
  public:

    using real_type            = realT;
    using parameter_type       = real_type;
    using container_type       = std::vector<parameter_type>;

    // `pair_parameter_type` is a parameter for a interacting pair.
    // Although it is the same type as `parameter_type` in this potential,
    // it can be different from normal parameter for each particle.
    // This enables NeighborList to cache a value to calculate forces between
    // the particles, e.g. by having (radius_i + radius_)/2 for a pair of {i, j}.
    using pair_parameter_type  = parameter_type;

    // ------------------------------------------------------------------------
    // topology stuff
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
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
        return parameter_type{0.0};
    }

  public:

    ExcludedVolumePotential(const real_type eps, const real_type cutoff_ratio,
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>&         exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : epsilon_(eps), cutoff_ratio_(cutoff_ratio),
          coef_at_cutoff_(std::pow(1.0 / cutoff_ratio, 12)),
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
    ~ExcludedVolumePotential() = default;
    ExcludedVolumePotential(const ExcludedVolumePotential&) = default;
    ExcludedVolumePotential(ExcludedVolumePotential&&)      = default;
    ExcludedVolumePotential& operator=(const ExcludedVolumePotential&) = default;
    ExcludedVolumePotential& operator=(ExcludedVolumePotential&&)      = default;

    // this value will be stored in NeighborList.
    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept
    {
        return this->parameters_[i] + this->parameters_[j];
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

    real_type potential(const real_type r, const pair_parameter_type& d) const noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
        MJOLNIR_LOG_FUNCTION_DEBUG();

        if(d * this->cutoff_ratio_ < r){return 0.0;}

        const real_type d_r  = d / r;
        const real_type dr3  = d_r * d_r * d_r;
        const real_type dr6  = dr3 * dr3;
        const real_type dr12 = dr6 * dr6;
        // correction by coef_at_cuttoff make energy 0 about the particle pair
        // who's distance is over the cutoff range.
        return this->epsilon_ * (dr12 - this->coef_at_cutoff_);
    }
    real_type derivative(const real_type r, const pair_parameter_type& d) const noexcept
    {
        if(d * this->cutoff_ratio_ < r){return 0.0;}

        const real_type rinv = 1.0 / r;
        const real_type d_r  = d * rinv;
        const real_type dr3  = d_r * d_r * d_r;
        const real_type dr6  = dr3 * dr3;
        const real_type dr12 = dr6 * dr6;
        return -12.0 * this->epsilon_ * dr12 * rinv;
    }

    template<typename traitsT>
    void initialize(const System<traitsT>& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->update(sys);
        return;
    }

    // nothing to be done if system parameter (e.g. temperature) changes
    template<typename traitsT>
    void update(const System<traitsT>& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // update exclusion list based on sys.topology()
        exclusion_list_.make(sys);
        return;
    }

    real_type cutoff_ratio()   const noexcept {return this->cutoff_ratio_;}
    real_type coef_at_cutoff() const noexcept {return this->coef_at_cutoff_;}

    real_type max_cutoff_length() const
    {
        const real_type max_sigma =
            *(std::max_element(parameters_.cbegin(), parameters_.cend()));
        return 2 * max_sigma * this->cutoff_ratio_;
    }

    // -----------------------------------------------------------------------
    // for spatial partitions

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
        return !exclusion_list_.is_excluded(i, j);
    }
    // for testing
    exclusion_list_type const& exclusion_list() const noexcept
    {
        return exclusion_list_;
    }

    // ------------------------------------------------------------------------
    // used by Observer.

    static const char* name() noexcept {return "ExcludedVolume";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    real_type& epsilon()       noexcept {return this->epsilon_;}
    real_type  epsilon() const noexcept {return this->epsilon_;}
    std::vector<real_type>&       parameters()       noexcept {return this->parameters_;}
    std::vector<real_type> const& parameters() const noexcept {return this->parameters_;}

  private:

    real_type epsilon_;
    real_type cutoff_ratio_;
    real_type coef_at_cutoff_; // correction of energy value at cutoff length
    std::vector<parameter_type> parameters_;
    std::vector<std::size_t> participants_;

    exclusion_list_type  exclusion_list_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class ExcludedVolumePotential<double>;
extern template class ExcludedVolumePotential<float>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif /* MJOLNIR_EXCLUDED_VOLUME_POTENTIAL */
