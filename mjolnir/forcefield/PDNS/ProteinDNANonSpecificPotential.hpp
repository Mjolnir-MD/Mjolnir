#ifndef MJOLNIR_FORCEFIELD_PDNS_PROTEIN_DNA_NON_SPECIFIC_POTENTIAL_HPP
#define MJOLNIR_FORCEFIELD_PDNS_PROTEIN_DNA_NON_SPECIFIC_POTENTIAL_HPP
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/core/System.hpp>
#include <cstdint>

namespace mjolnir
{

// Protein-DNA Non-Specific interaction that represents hydrogen bond between
// side chain atoms in proteins and the phosphate residues in DNA.
// This is an implementation of the potential developed in the following paper.
// - T.Niina, G.B.Brandani, C.Tan, and S.Takada (2017) PLoS. Comput. Biol.
//
// ```toml
// [[forcefields.global]]
// interaction = "PDNS"
// potential   = "PDNS"
// spatial_partition.type   = "CellList"
// spatial_partition.margin = 1.0
// sigma = 1.0
// delta = 0.17453
// parameters  = [
// {index =    2, S3 = 1, kind = "DNA"},
// {index =    5, S3 = 4, kind = "DNA"},
// # ...
// {index = 1000, kind = "Protein", PN =  999, PC = 1001, k = 1.2, r0 = 5.0, theta0 = 100.0, phi0 = 130.0},
// {index = 1023, kind = "Protein", PN = 1022, PC = 1024, k = 1.2, r0 = 6.0, theta0 = 110.0, phi0 = 120.0},
// # ...
// ]
// ```

template<typename realT>
class ProteinDNANonSpecificPotential
{
  public:
    using real_type = realT;
    using self_type = ProteinDNANonSpecificPotential<real_type>;

    struct parameter_type
    {
        std::uint32_t P, PN, PC; // for PRO
        real_type k, r0, theta0, phi0;
        real_type r_cut, r_cut_sq;
    };
    using container_type = std::vector<parameter_type>;
    using dna_index_type = std::pair<std::uint32_t, std::uint32_t>;

    static constexpr std::uint32_t invalid() noexcept
    {
        return std::numeric_limits<std::uint32_t>::max();
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{invalid(), invalid(), 0, 0, 0, 0};
    }
    static constexpr real_type default_cutoff() noexcept {return 5.0;}

  public:

    ProteinDNANonSpecificPotential(const real_type sigma,
        const real_type delta, const real_type cutoff_ratio,
        const std::vector<parameter_type>& parameters,
        const std::vector<dna_index_type>& dna_idxs)
        : sigma_(sigma), delta_(delta), delta2_(delta * 2),
          pi_over_2delta_(math::constants<real_type>::pi() * 0.5 / delta),
          cutoff_ratio_(cutoff_ratio), max_cutoff_length_(0),
          parameters_(parameters), dnas_(dna_idxs)
    {}
    ~ProteinDNANonSpecificPotential() = default;

    template<typename traitsT>
    void initialize(const System<traitsT>& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->update(sys);
        return;
    }

    template<typename traitsT>
    void update(const System<traitsT>&) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        // re-set cutoff length

        this->max_cutoff_length_ = real_type(0);
        for(auto& para : this->parameters_)
        {
            para.r_cut    = para.r0 + cutoff_ratio_ * sigma_;
            para.r_cut_sq = para.r_cut * para.r_cut;

            this->max_cutoff_length_ = std::max(max_cutoff_length_, para.r_cut);
        }
        return;
    }

    std::pair<real_type, real_type>
    f_df(const real_type r0, const real_type r) const noexcept
    {
        const real_type rsigma   = real_type(1) / this->sigma_;
        const real_type dr_sigma = (r - r0) * rsigma;
        const real_type term     = std::exp(real_type(-0.5) * dr_sigma * dr_sigma);

        return std::make_pair(term, -dr_sigma * rsigma * term);
    }
    std::pair<real_type, real_type>
    g_dg(const real_type theta0, const real_type theta) const noexcept
    {
        const real_type dtheta     = theta - theta0;
        const real_type abs_dtheta = std::abs(dtheta);
        if(abs_dtheta < this->delta_)
        {
            return std::make_pair(real_type(1), real_type(0));
        }
        else if(abs_dtheta < this->delta2_)
        {
            const real_type c = std::cos(dtheta * pi_over_2delta_);
            const real_type s = std::sin(dtheta * pi_over_2delta_);
            return std::make_pair(1 - c * c, 2 * s * c * pi_over_2delta_);
        }
        else
        {
            return std::make_pair(real_type(0), real_type(0));
        }
    }

    real_type f(const real_type r0, const real_type r) const noexcept
    {
        const real_type dr_sigma = (r - r0) / this->sigma_;
        return std::exp(real_type(-0.5) * dr_sigma * dr_sigma);
    }
    real_type g(const real_type theta0, const real_type theta) const noexcept
    {
        const real_type dtheta     = theta - theta0;
        const real_type abs_dtheta = std::abs(dtheta);
        if(abs_dtheta < this->delta_)
        {
            return 1;
        }
        else if(abs_dtheta < this->delta2_)
        {
            const real_type term = std::cos(dtheta * pi_over_2delta_);
            return real_type(1) - term * term;
        }
        else
        {
            return 0;
        }
    }

    real_type cutoff_ratio()      const noexcept {return cutoff_ratio_;}
    real_type max_cutoff_length() const noexcept {return this->max_cutoff_length_;}

    // access to the parameters...
    std::vector<parameter_type>&       contacts()       noexcept {return parameters_;}
    std::vector<parameter_type> const& contacts() const noexcept {return parameters_;}
    std::vector<dna_index_type>&       dnas()       noexcept {return dnas_;}
    std::vector<dna_index_type> const& dnas() const noexcept {return dnas_;}

    // ------------------------------------------------------------------------
    // used by Observer.
    static const char* name() noexcept {return "PDNS";}

  private:

    real_type sigma_, delta_, delta2_, pi_over_2delta_;
    real_type cutoff_ratio_, max_cutoff_length_;
    std::vector<parameter_type> parameters_;
    std::vector<dna_index_type> dnas_;
};

} // mjolnir
#endif// MJOLNIR_FORCEFIELD_PDNS_PROTEIN_DNA_NON_SPECIFIC_INTERACTION_HPP
