#ifndef MJOLNIR_EXCLUDED_VOLUME_POTENTIAL
#define MJOLNIR_EXCLUDED_VOLUME_POTENTIAL
#include <mjolnir/core/System.hpp>
#include <algorithm>
#include <cmath>

namespace mjolnir
{

/*! @brief excluded volume potential        *
 *  V(r) = epsilon * (sigma/r)^12           *
 * dV/dr = -12 * epsilon * (sigma/r)^12 / r */
template<typename traitsT>
class ExcludedVolumePotential
{
  public:

    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef real_type parameter_type;
    typedef std::vector<parameter_type> container_type;

    // topology stuff
    typedef StructureTopology topology_type;
    typedef typename topology_type::group_id_type        group_id_type;
    typedef typename topology_type::connection_name_type connection_name_type;

    // rc = 2.0 * sigma
    constexpr static real_type cutoff_ratio = 2.0;

  public:

    ExcludedVolumePotential(const real_type eps, const container_type& params)
        : epsilon_(eps), radii_(params)
    {}
    ExcludedVolumePotential(const real_type eps, container_type&& params)
        : epsilon_(eps), radii_(std::move(params))
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
        return this->epsilon_ * dr12;
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
    void update(const System<traitsT>&) const noexcept {return;}

    real_type max_cutoff_length() const
    {
        const real_type max_sigma =
            *(std::max_element(radii_.cbegin(), radii_.cend()));
        return 2 * max_sigma * cutoff_ratio;
    }

    // e.g. {"bond", 3} means ignore particles connected within 3 "bond"s
    std::vector<std::pair<connection_name_type, std::size_t>>&
    ignored_connections()       noexcept {return this->ignored_connections_;};
    std::vector<std::pair<connection_name_type, std::size_t>> const&
    ignored_connections() const noexcept {return this->ignored_connections_;};

    bool is_ignored_group(const group_id_type& i, const group_id_type& j) const
    {
        if     (this->ignored_group_.empty())     {return false;}
        else if(this->ignored_group_ == "self")   {return i == j;}
        else if(this->ignored_group_ == "others") {return i != j;}

        throw std::runtime_error(
                "unknown group ignoration option: " + this->ignored_group_);
    }

    std::string name() const noexcept {return "ExcludedVolume";}

    // access to the parameters
    real_type& epsilon()       noexcept {return this->epsilon_;}
    real_type  epsilon() const noexcept {return this->epsilon_;}
    std::vector<real_type>&       radii()       noexcept {return this->radii_;}
    std::vector<real_type> const& radii() const noexcept {return this->radii_;}

  private:

    real_type epsilon_;
    std::vector<parameter_type> radii_;

    std::string ignored_group_;
    std::vector<std::pair<connection_name_type, std::size_t>> ignored_connections_;
};

template<typename traitsT>
constexpr typename ExcludedVolumePotential<traitsT>::real_type
ExcludedVolumePotential<traitsT>::cutoff_ratio;

} // mjolnir
#endif /* MJOLNIR_EXCLUDED_VOLUME_POTENTIAL */
