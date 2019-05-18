#ifndef MJOLNIR_POTENTIAL_GLOBAL_EXCLUDED_VOLUME_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_EXCLUDED_VOLUME_POTENTIAL_HPP
#include <mjolnir/potential/global/IgnoreMolecule.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/math.hpp>
#include <algorithm>
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
template<typename realT>
class ExcludedVolumePotential
{
  public:

    using real_type      = realT;
    using parameter_type = real_type;
    using container_type = std::vector<parameter_type>;

    // topology stuff
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;

    // rc = 2.0 * sigma
    constexpr static real_type cutoff_ratio = 2.0;

    // to make the potential curve continuous at the cutoff point
    constexpr static real_type coef_at_cutoff =
        compiletime::pow(1.0 / cutoff_ratio, 12);

  public:

    ExcludedVolumePotential(const real_type eps, container_type params,
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_molecule)
        : epsilon_(eps),
          radii_(std::move(params)),
          ignore_molecule_(std::move(ignore_molecule)),
          ignore_within_  (exclusions.begin(), exclusions.end())
    {}
    ~ExcludedVolumePotential() = default;
    ExcludedVolumePotential(const ExcludedVolumePotential&) = default;
    ExcludedVolumePotential(ExcludedVolumePotential&&)      = default;
    ExcludedVolumePotential& operator=(const ExcludedVolumePotential&) = default;
    ExcludedVolumePotential& operator=(ExcludedVolumePotential&&)      = default;

    parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept
    {
        return this->radii_[i] + this->radii_[j];
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

    real_type potential(const real_type r, const parameter_type& d) const noexcept
    {
        if(d * cutoff_ratio < r){return 0.0;}

        const real_type d_r  = d / r;
        const real_type dr3  = d_r * d_r * d_r;
        const real_type dr6  = dr3 * dr3;
        const real_type dr12 = dr6 * dr6;
        return this->epsilon_ * (dr12 - coef_at_cutoff);
    }
    real_type derivative(const real_type r, const parameter_type& d) const noexcept
    {
        if(d * cutoff_ratio < r){return 0.0;}

        const real_type rinv = 1.0 / r;
        const real_type d_r  = d * rinv;
        const real_type dr3  = d_r * d_r * d_r;
        const real_type dr6  = dr3 * dr3;
        const real_type dr12 = dr6 * dr6;
        return -12.0 * this->epsilon_ * dr12 * rinv;
    }

    template<typename traitsT>
    void initialize(const System<traitsT>&) const noexcept {return;}

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

    bool is_ignored_molecule(
            const molecule_id_type& i, const molecule_id_type& j) const noexcept
    {
        return ignore_molecule_.is_ignored(i, j);
    }

    static const char* name() noexcept {return "ExcludedVolume";}

    // access to the parameters
    real_type& epsilon()       noexcept {return this->epsilon_;}
    real_type  epsilon() const noexcept {return this->epsilon_;}
    std::vector<real_type>&       parameters()       noexcept {return this->radii_;}
    std::vector<real_type> const& parameters() const noexcept {return this->radii_;}

  private:

    real_type epsilon_;
    std::vector<real_type> radii_;

    ignore_molecule_type ignore_molecule_;
    std::vector<std::pair<connection_kind_type, std::size_t>> ignore_within_;
};

template<typename realT>
constexpr typename ExcludedVolumePotential<realT>::real_type
ExcludedVolumePotential<realT>::cutoff_ratio;

} // mjolnir
#endif /* MJOLNIR_EXCLUDED_VOLUME_POTENTIAL */
