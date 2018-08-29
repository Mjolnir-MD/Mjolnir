#ifndef MJOLNIR_LENNARD_JONES_POTENTIAL
#define MJOLNIR_LENNARD_JONES_POTENTIAL
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/potential/global/ChainIgnoration.hpp>
#include <vector>
#include <algorithm>
#include <cmath>

namespace mjolnir
{

/*! @brief Lennard-Jones type potential & derivative                       *
 * designed for global force field.                                        *
 * V(r)  =  4. * epsilon * ((r/sigma)^12 - (r/sigma)^6))                   *
 * dV/dr = 24. * epsilon / r * ((r/sigma)^6 - 2 * (r/sigma)^12)            */
template<typename realT, typename ChainIgnoration>
class LennardJonesPotential
{
  public:
    using real_type = realT;
    using parameter_type = std::pair<real_type, real_type>; // {sigma, epsilon}
    using container_type = std::vector<parameter_type>;
    // topology stuff
    using topology_type  = Topology;
    using chain_id_type  = typename topology_type::chain_id_type;
    using connection_kind_type  = typename topology_type::connection_kind_type;
    using chain_ignoration_type = ChainIgnoration;

    // rc = 2.5 * sigma
    constexpr static real_type cutoff_ratio = 2.5;
    // to make the potential curve continuous at the cutoff point
    constexpr static real_type coef_at_cutoff =
        compiletime::pow(1.0 / cutoff_ratio, 12u) -
        compiletime::pow(1.0 / cutoff_ratio,  6u);

  public:

    LennardJonesPotential(const std::vector<parameter_type>& radii,
        const std::map<connection_kind_type, std::size_t>& exclusions)
        : radii_(radii), ignored_chain_(),
          ignore_within_(exclusions.begin(), exclusions.end())
    {}
    LennardJonesPotential(std::vector<parameter_type>&& radii,
        const std::map<connection_kind_type, std::size_t>& exclusions)
        : radii_(std::move(radii)), ignored_chain_(),
          ignore_within_(exclusions.begin(), exclusions.end())
    {}
    ~LennardJonesPotential() = default;

    real_type potential(const std::size_t i, const std::size_t j,
                        const real_type r) const noexcept
    {
        const real_type sigma = 0.5 * (radii_[i].first + radii_[j].first);
        if(sigma * cutoff_ratio < r){return 0.0;}

        const real_type epsilon = (radii_[i].second == radii_[j].second) ?
           radii_[i].second : std::sqrt(radii_[i].second * radii_[j].second);

        const real_type r1s1   = sigma / r;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return 4.0 * epsilon * (r12s12 - r6s6 - coef_at_cutoff);
    }
    real_type potential(const std::size_t i, const std::size_t j,
                        const real_type   r, const parameter_type& p) const noexcept
    {
        const real_type sigma = p.first;
        if(sigma * cutoff_ratio < r){return 0.0;}

        const real_type epsilon = p.second;

        const real_type r1s1   = sigma / r;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return 4.0 * epsilon * (r12s12 - r6s6 - coef_at_cutoff);
    }

    real_type derivative(const std::size_t i, const std::size_t j,
                         const real_type r) const noexcept
    {
        const real_type sigma = 0.5 * (radii_[i].first + radii_[j].first);
        if(sigma * cutoff_ratio < r){return 0.0;}

        const real_type epsilon = (radii_[i].second == radii_[j].second) ?
           radii_[i].second : std::sqrt(radii_[i].second * radii_[j].second);

        const real_type r1s1   = sigma / r;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return 24.0 * epsilon * (r6s6 - 2.0 * r12s12) / r;
    }
    real_type derivative(const std::size_t i, const std::size_t j,
                         const real_type   r, const parameter_type& p) const noexcept
    {
        const real_type sigma = p.first;
        if(sigma * cutoff_ratio < r){return 0.0;}

        const real_type epsilon = p.second;

        const real_type r1s1   = sigma / r;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return 24.0 * epsilon * (r6s6 - 2.0 * r12s12) / r;
    }

    parameter_type prepair_params(std::size_t i, std::size_t j) const noexcept
    {
        const auto sgm1 = radii_[i].first;
        const auto eps1 = radii_[i].second;
        const auto sgm2 = radii_[j].first;
        const auto eps2 = radii_[j].second;

        return std::make_pair((sgm1 + sgm2) / 2,
                             ((eps1 == eps2) ? eps1 : std::sqrt(eps1 * eps2)));
    }

    real_type max_cutoff_length() const noexcept
    {
        const real_type max_sigma = std::max_element(
            this->radii_.cbegin(), this->radii_.cend(),
            [](const parameter_type& lhs, const parameter_type& rhs) noexcept {
                return lhs.first < rhs.first;
            })->first;
        return max_sigma * cutoff_ratio;
    }

    // nothing to do when system parameters change.
    template<typename traitsT>
    void update(const System<traitsT>& sys) const noexcept {return;}

    // e.g. `3` means ignore particles connected within 3 "bond"s
    std::vector<std::pair<connection_kind_type, std::size_t>>
    ignore_within() const {return ignore_within_;}

    bool is_ignored_chain(
            const chain_id_type& i, const chain_id_type& j) const noexcept
    {
        return ignored_chain_.is_ignored(i, j);
    }

    static const char* name() noexcept {return "LennardJones";}

    // access to the parameters...
    std::vector<parameter_type>&       radii()       noexcept {return radii_;}
    std::vector<parameter_type> const& radii() const noexcept {return radii_;}

  private:

    container_type radii_;

    chain_ignoration_type ignored_chain_;
    std::vector<std::pair<connection_kind_type, std::size_t>> ignore_within_;
};
template<typename realT, typename ignoreT>
constexpr typename LennardJonesPotential<realT, ignoreT>::real_type
LennardJonesPotential<realT, ignoreT>::cutoff_ratio;

} // mjolnir
#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
