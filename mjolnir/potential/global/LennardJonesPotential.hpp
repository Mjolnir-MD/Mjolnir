#ifndef MJOLNIR_POTENTIAL_GLOBAL_LENNARD_JONES_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_LENNARD_JONES_POTENTIAL_HPP
#include <mjolnir/potential/global/IgnoreMolecule.hpp>
#include <mjolnir/potential/global/IgnoreGroup.hpp>
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
    using real_type = realT;
    using parameter_type = std::pair<real_type, real_type>; // {sigma, epsilon}
    using container_type = std::vector<parameter_type>;

    // topology stuff
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;

    // rc = 2.5 * sigma
    constexpr static real_type cutoff_ratio = 2.5;
    // to make the potential curve continuous at the cutoff point
    constexpr static real_type coef_at_cutoff =
        compiletime::pow<real_type>(1 / cutoff_ratio, 12u) -
        compiletime::pow<real_type>(1 / cutoff_ratio,  6u);

    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{real_type(0), real_type(0)};
    }

  public:

    LennardJonesPotential(
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>&         exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : ignore_molecule_(std::move(ignore_mol)),
          ignore_group_   (std::move(ignore_grp)),
          ignore_within_  (exclusions.begin(), exclusions.end())
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
    ~LennardJonesPotential() = default;

    parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept
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

    real_type potential(const real_type r, const parameter_type& p) const noexcept
    {
        const real_type sigma = p.first;
        if(sigma * cutoff_ratio < r){return 0;}

        const real_type epsilon = p.second;

        const real_type r1s1   = sigma / r;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return 4 * epsilon * (r12s12 - r6s6 - coef_at_cutoff);
    }
    real_type derivative(const real_type r, const parameter_type& p) const noexcept
    {
        const real_type sigma = p.first;
        if(sigma * cutoff_ratio < r){return 0;}

        const real_type epsilon = p.second;

        const real_type rinv   = 1 / r;
        const real_type r1s1   = sigma * rinv;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return 24 * epsilon * (r6s6 - 2 * r12s12) * rinv;
    }

    real_type max_cutoff_length() const noexcept
    {
        const real_type max_sigma = std::max_element(
            this->parameters_.cbegin(), this->parameters_.cend(),
            [](const parameter_type& lhs, const parameter_type& rhs) noexcept {
                return lhs.first < rhs.first;
            })->first;
        return max_sigma * cutoff_ratio;
    }

    template<typename traitsT>
    void initialize(const System<traitsT>&) noexcept {return;}

    // nothing to do when system parameters change.
    template<typename traitsT>
    void update(const System<traitsT>&) const noexcept {return;}

    // ------------------------------------------------------------------------
    // ignore_xxx functions would be used in core/ExclusionLists.

    // e.g. `{"bond", 3}` means ignore particles connected within 3 "bond"s
    std::vector<std::pair<connection_kind_type, std::size_t>>
    ignore_within() const {return ignore_within_;}

    bool is_ignored_molecule(
            const molecule_id_type& i, const molecule_id_type& j) const noexcept
    {
        return ignore_molecule_.is_ignored(i, j);
    }
    bool is_ignored_group(
            const group_id_type& i, const group_id_type& j) const noexcept
    {
        return ignore_group_.is_ignored(i, j);
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

    container_type parameters_;
    std::vector<std::size_t> participants_;

    ignore_molecule_type ignore_molecule_;
    ignore_group_type    ignore_group_;
    std::vector<std::pair<connection_kind_type, std::size_t>> ignore_within_;
};
template<typename realT>
constexpr typename LennardJonesPotential<realT>::real_type
LennardJonesPotential<realT>::cutoff_ratio;

} // mjolnir
#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
