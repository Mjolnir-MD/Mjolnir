#ifndef MJOLNIR_POTENTIAL_GLOBAL_LENNARD_JONES_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_LENNARD_JONES_POTENTIAL_HPP
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/math.hpp>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace mjolnir
{

// Well-known Lennard-Jones interaction with Lorentz-Berthelot combining rules.
// This class contains sigmas and epsilons of the particles and calculates
// energy and derivative of the potential function.
template<typename realT>
class LennardJonesPotential
{
  public:
    using real_type            = realT;
    using parameter_type       = std::pair<real_type, real_type>; // {sigma, epsilon}
    using container_type       = std::vector<parameter_type>;

    // `pair_parameter_type` is a parameter for a interacting pair.
    // Although it is the same type as `parameter_type` in this potential,
    // it can be different from normal parameter for each particle.
    // This enables NeighborList to cache a value to calculate forces between
    // the particles, e.g. by having (sigma_i + sigma_j)/2 for a pair of {i, j}.
    using pair_parameter_type  = parameter_type;

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
        return real_type(2.5);
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{real_type(0), real_type(0)};
    }

  public:

    LennardJonesPotential(const real_type cutoff_ratio,
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>&         exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
    : cutoff_ratio_(cutoff_ratio),
      coef_at_cutoff_(std::pow(1 / cutoff_ratio, 12) - std::pow(1 / cutoff_ratio, 6)),
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
    LennardJonesPotential(
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>&         exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
    : LennardJonesPotential(default_cutoff(), parameters, exclusions,
                            std::move(ignore_mol), std::move(ignore_grp))
    {}
    ~LennardJonesPotential() = default;

    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept
    {
        const auto sgm1 = parameters_[i].first;
        const auto eps1 = parameters_[i].second;
        const auto sgm2 = parameters_[j].first;
        const auto eps2 = parameters_[j].second;

        return std::make_pair((sgm1 + sgm2) / 2,
                             ((eps1 == eps2) ? eps1 : std::sqrt(eps1 * eps2)));
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
        if(sigma * this->cutoff_ratio_ < r){return 0;}

        const real_type epsilon = p.second;

        const real_type r1s1   = sigma / r;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return 4 * epsilon * (r12s12 - r6s6 - coef_at_cutoff_);
    }
    real_type derivative(const real_type r, const pair_parameter_type& p) const noexcept
    {
        const real_type sigma = p.first;
        if(sigma * this->cutoff_ratio_ < r){return 0;}

        const real_type epsilon = p.second;

        const real_type rinv   = 1 / r;
        const real_type r1s1   = sigma * rinv;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return 24 * epsilon * (r6s6 - 2 * r12s12) * rinv;
    }

    real_type cutoff_ratio()   const noexcept {return this->cutoff_ratio_;}
    real_type coef_at_cutoff() const noexcept {return this->coef_at_cutoff_;}

    real_type max_cutoff_length() const noexcept
    {
        const real_type max_sigma = std::max_element(
            this->parameters_.cbegin(), this->parameters_.cend(),
            [](const parameter_type& lhs, const parameter_type& rhs) noexcept {
                return lhs.first < rhs.first;
            })->first;
        return max_sigma * this->cutoff_ratio_;
    }

    template<typename traitsT>
    void initialize(const System<traitsT>& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->update(sys);
        return;
    }

    template<typename traitsT>
    void update(const System<traitsT>& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // update exclusion list based on sys.topology()
        exclusion_list_.make(sys);
        return;
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
    static const char* name() noexcept {return "LennardJones";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    // access to the parameters...
    std::vector<parameter_type>&       parameters()       noexcept {return parameters_;}
    std::vector<parameter_type> const& parameters() const noexcept {return parameters_;}

    std::vector<std::size_t> const& participants() const noexcept {return participants_;}

  private:

    real_type cutoff_ratio_;
    real_type coef_at_cutoff_;
    container_type parameters_;
    std::vector<std::size_t> participants_;

    exclusion_list_type  exclusion_list_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class LennardJonesPotential<double>;
extern template class LennardJonesPotential<float>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
