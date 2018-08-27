#ifndef MJOLNIR_EXCLUDED_VOLUME_POTENTIAL
#define MJOLNIR_EXCLUDED_VOLUME_POTENTIAL
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/potential/global/ChainIgnoration.hpp>
#include <algorithm>
#include <cmath>

namespace mjolnir
{

/*! @brief excluded volume potential        *
 *  V(r) = epsilon * (sigma/r)^12           *
 * dV/dr = -12 * epsilon * (sigma/r)^12 / r */
template<typename realT, typename ChainIgnoration>
class ExcludedVolumePotential
{
  public:

    using real_type      = realT;
    using parameter_type = real_type;
    using container_type = std::vector<parameter_type>;

    // topology stuff
    using topology_type = Topology;
    using chain_id_type = typename topology_type::chain_id_type;
    using connection_kind_type  = typename topology_type::connection_kind_type;
    using chain_ignoration_type = ChainIgnoration;

    // rc = 2.0 * sigma
    constexpr static real_type cutoff_ratio = 2.0;

    // to make the potential curve continuous at the cutoff point
    constexpr static real_type coef_at_cutoff =
        compiletime::pow(1.0 / cutoff_ratio, 12);

  public:

    ExcludedVolumePotential(const real_type eps, const container_type& params,
        const std::map<connection_kind_type, std::size_t>& exclusions)
        : epsilon_(eps), radii_(params), ignored_chain_(),
          ignore_within_(exclusions.begin(), exclusions.end())
    {}
    ExcludedVolumePotential(const real_type eps, container_type&& params,
        const std::map<connection_kind_type, std::size_t>& exclusions)
        : epsilon_(eps), radii_(std::move(params)), ignored_chain_(),
          ignore_within_(exclusions.begin(), exclusions.end())
    {}
    ~ExcludedVolumePotential() = default;

    ExcludedVolumePotential(const ExcludedVolumePotential&) = default;
    ExcludedVolumePotential(ExcludedVolumePotential&&)      = default;
    ExcludedVolumePotential& operator=(const ExcludedVolumePotential&) = default;
    ExcludedVolumePotential& operator=(ExcludedVolumePotential&&)      = default;

    real_type potential (const std::size_t i, const std::size_t j,
                         const real_type r) const noexcept
    {
        const real_type d = this->radii_[i] + this->radii_[j];
        if(d * cutoff_ratio < r){return 0.0;}

        const real_type d_r  = d / r;
        const real_type dr3  = d_r * d_r * d_r;
        const real_type dr6  = dr3 * dr3;
        const real_type dr12 = dr6 * dr6;
        return this->epsilon_ * (dr12 - coef_at_cutoff);
    }

    real_type derivative(const std::size_t i, const std::size_t j,
                         const real_type r) const noexcept
    {
        const real_type d = this->radii_[i] + this->radii_[j];
        if(d * cutoff_ratio < r){return 0.0;}

        const real_type rinv = 1. / r;
        const real_type d_r  = d * rinv;
        const real_type dr3  = d_r * d_r * d_r;
        const real_type dr6  = dr3 * dr3;
        const real_type dr12 = dr6 * dr6;
        return -12.0 * this->epsilon_ * dr12 * rinv;
    }

    // nothing to be done if system parameter (e.g. temperature) changes
    template<typename traitsT>
    void update(const System<traitsT>&) const noexcept {return;}

    real_type max_cutoff_length() const
    {
        const real_type max_sigma =
            *(std::max_element(radii_.cbegin(), radii_.cend()));
        return 2 * max_sigma * cutoff_ratio;
    }

    // e.g. "bond" -> 3 means ignore particles connected within 3 "bond"s
    std::vector<std::pair<connection_kind_type, std::size_t>>
    ignore_within() const {return ignore_within_;}

    bool is_ignored_chain(
            const chain_id_type& i, const chain_id_type& j) const noexcept
    {
        return ignored_chain_.is_ignored(i, j);
    }

    static const char* name() noexcept {return "ExcludedVolume";}

    // access to the parameters
    real_type& epsilon()       noexcept {return this->epsilon_;}
    real_type  epsilon() const noexcept {return this->epsilon_;}
    std::vector<real_type>&       radii()       noexcept {return this->radii_;}
    std::vector<real_type> const& radii() const noexcept {return this->radii_;}

  private:

    real_type epsilon_;
    std::vector<real_type> radii_;

    chain_ignoration_type ignored_chain_;
    std::vector<std::pair<connection_kind_type, std::size_t>> ignore_within_;
};

template<typename realT, typename ignoreT>
constexpr typename ExcludedVolumePotential<realT, ignoreT>::real_type
ExcludedVolumePotential<realT, ignoreT>::cutoff_ratio;

} // mjolnir
#endif /* MJOLNIR_EXCLUDED_VOLUME_POTENTIAL */
