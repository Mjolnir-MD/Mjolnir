#ifndef MJOLNIR_LENNARD_JONES_POTENTIAL
#define MJOLNIR_LENNARD_JONES_POTENTIAL
#include <mjolnir/core/System.hpp>
#include <mjolnir/potential/global/ChainIgnoration.hpp>
#include <vector>
#include <algorithm>
#include <cmath>

namespace mjolnir
{

/*! @brief Lennard-Jones type potential & derivative                       *
 * designed for global force field. so it doesn't have its own parameters. *
 * V(r)  =  4. * epsilon * ((r/sigma)^12 - (r/sigma)^6))                   *
 * dV/dr = 24. * epsilon / r * ((r/sigma)^6 - 2 * (r/sigma)^12)            */
template<typename traitsT, typename ChainIgnoration>
class LennardJonesPotential
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef std::pair<real_type, real_type> parameter_type;

    // topology stuff
    typedef StructureTopology topology_type;
    typedef typename topology_type::chain_id_type        chain_id_type;
    typedef typename topology_type::connection_kind_type connection_kind_type;
    typedef ChainIgnoration chain_ignoration_type;

    // rc = 2.5 * sigma
    constexpr static real_type cutoff_ratio = 2.5;

  public:
    LennardJonesPotential() = default;
    LennardJonesPotential(const std::vector<parameter_type>& radii)
        : radii_(radii)
    {}
    LennardJonesPotential(std::vector<parameter_type>&& radii)
        : radii_(std::forward<std::vector<parameter_type>>(radii))
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
        return 4.0 * epsilon * (r12s12 - r6s6);
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
    void update(const System<traitsT>& sys) const noexcept {return;}

    // e.g. `3` means ignore particles connected within 3 "bond"s
    std::size_t  ignored_bonds()    const noexcept {return ignored_bonds_;}
    std::size_t& ignored_bonds()          noexcept {return ignored_bonds_;}
    std::size_t  ignored_contacts() const noexcept {return ignored_contacts_;}
    std::size_t& ignored_contacts()       noexcept {return ignored_contacts_;}

    bool is_ignored_chain(
            const chain_id_type& i, const chain_id_type& j) const noexcept
    {
        return ignored_chain_.is_ignored(i, j);
    }

    std::string name() const noexcept {return "LennardJones";}

    // access to the parameters...
    std::vector<parameter_type>&       radii()       noexcept {return radii_;}
    std::vector<parameter_type> const& radii() const noexcept {return radii_;}

  private:

    std::vector<parameter_type> radii_;

    chain_ignoration_type ignored_chain_;
    std::size_t ignored_bonds_;
    std::size_t ignored_contacts_;
};
template<typename traitsT, typename ignoreT>
constexpr typename LennardJonesPotential<traitsT, ignoreT>::real_type
LennardJonesPotential<traitsT, ignoreT>::cutoff_ratio;

} // mjolnir
#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
