#ifndef MJOLNIR_FORCEFIELD_PDNS_PROTEIN_DNA_NON_SPECIFIC_POTENTIAL_HPP
#define MJOLNIR_FORCEFIELD_PDNS_PROTEIN_DNA_NON_SPECIFIC_POTENTIAL_HPP
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/Topology.hpp>
#include <cstdint>

namespace mjolnir
{

// Protein-DNA Non-Specific interaction that represents hydrogen bond between
// side chain atoms in proteins and the phosphate residues in DNA.
// This is an implementation of the potential developed in the following paper.
// - T.Niina, G.B.Brandani, C.Tan, and S.Takada (2017) PLoS. Comput. Biol.

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
namespace parameter_PDNS
{
enum bead_kind : std::uint8_t
{
    Protein,
    DNA,
    invalid
};

template<typename charT, typename traits>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const bead_kind bk)
{
    switch(bk)
    {
        case bead_kind::Protein: {os << "PRO"; return os;}
        case bead_kind::DNA:     {os << "DNA"; return os;}
        case bead_kind::invalid: {os << "INV"; return os;}
        default:                 {os << "UNK"; return os;}
    }
}
} // pdns

template<typename realT>
class ProteinDNANonSpecificPotential
{
  public:
    using real_type = realT;
    using self_type = ProteinDNANonSpecificPotential<real_type>;
    using bead_kind = parameter_PDNS::bead_kind;

    struct parameter_type
    {
        bead_kind     kind;
        std::uint32_t S3;     // for DNA
        std::uint32_t PN, PC; // for PRO
        real_type k, r0, theta0, phi0;
    };
    struct pair_parameter_type
    {
        std::uint32_t S3, PN, PC, DNA; // `DNA` represents which idx is DNA
        real_type k, r0, theta0, phi0;
    };
    using container_type = std::vector<parameter_type>;

    // topology stuff
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList;

    static constexpr std::uint32_t invalid() noexcept
    {
        return std::numeric_limits<std::uint32_t>::max();
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{
            bead_kind::invalid, invalid(), invalid(), invalid(), 0, 0, 0, 0
        };
    }
    static constexpr real_type default_cutoff() noexcept {return 5.0;}

  public:

    ProteinDNANonSpecificPotential(const real_type sigma, const real_type delta,
        const real_type cutoff_ratio,
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>&         exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : sigma_(sigma), delta_(delta), delta2_(delta * 2),
          one_over_2delta_(0.5 / delta), cutoff_ratio_(cutoff_ratio),
          exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->parameters_  .reserve(parameters.size());
        this->participants_.reserve(parameters.size());
        for(const auto& idxp : parameters)
        {
            const auto idx = idxp.first;

            this->participants_.push_back(idx);
            if(idxp.second.kind == bead_kind::Protein)
            {
                this->proteins_.push_back(idx);
            }

            if(idx >= this->parameters_.size())
            {
                this->parameters_.resize(idx+1, self_type::default_parameter());
            }
            this->parameters_.at(idx) = idxp.second;
        }
    }
    ~ProteinDNANonSpecificPotential() = default;

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
        constexpr auto  pi         = math::constants<real_type>::pi();
        const real_type dtheta     = theta - theta0;
        const real_type abs_dtheta = std::abs(dtheta);
        if(abs_dtheta < this->delta_)
        {
            return std::make_pair(real_type(1), real_type(0));
        }
        else if(abs_dtheta < this->delta2_)
        {
            const real_type c = std::cos(pi * dtheta * one_over_2delta_);
            const real_type s = std::sin(pi * dtheta * one_over_2delta_);
            return std::make_pair(1 - c * c, 2 * s * c * pi * one_over_2delta_);
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
        constexpr auto pi = math::constants<real_type>::pi();
        const real_type dtheta     = theta - theta0;
        const real_type abs_dtheta = std::abs(dtheta);
        if(abs_dtheta < this->delta_)
        {
            return 1;
        }
        else if(abs_dtheta < this->delta2_)
        {
            const real_type term = std::cos(pi * dtheta * one_over_2delta_);
            return real_type(1) - term * term;
        }
        else
        {
            return 0;
        }
    }

    pair_parameter_type
    prepare_params(const std::size_t i, const std::size_t j) const noexcept
    {
        const auto& pi = this->parameters_[i];
        const auto& pj = this->parameters_[j];
        assert(pi.kind != bead_kind::invalid);
        assert(pj.kind != bead_kind::invalid);
        assert(pi.kind != pj.kind);

        pair_parameter_type pp;

        if(pi.kind == bead_kind::DNA && pj.kind == bead_kind::Protein)
        {
            pp.DNA    = static_cast<std::uint32_t>(i);
            pp.S3     = pi.S3;
            pp.PN     = pj.PN;
            pp.PC     = pj.PC;
            pp.k      = pj.k;
            pp.r0     = pj.r0;
            pp.theta0 = pj.theta0;
            pp.phi0   = pj.phi0;
        }
        else if(pi.kind == bead_kind::Protein && pj.kind == bead_kind::DNA)
        {
            pp.DNA    = static_cast<std::uint32_t>(j);
            pp.S3     = pj.S3;
            pp.PN     = pi.PN;
            pp.PC     = pi.PC;
            pp.k      = pi.k;
            pp.r0     = pi.r0;
            pp.theta0 = pi.theta0;
            pp.phi0   = pi.phi0;
        }
        else
        {
            assert(false);
        }
        return pp;
    }

    real_type cutoff_ratio()      const noexcept {return cutoff_ratio_;}
    real_type gaussian_cutoff()   const noexcept {return sigma_ * cutoff_ratio_;}
    real_type max_cutoff_length() const noexcept
    {
        real_type max_r0 = 0.0;
        for(const auto pro : this->proteins_)
        {
            const auto& para = this->parameters_.at(pro);
            max_r0 = std::max(max_r0, para.r0);
        }
        return max_r0 + this->sigma_ * this->cutoff_ratio_;
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
        if(exclusion_list_.is_excluded(i, j))
        {
            return false;
        }
        const auto i_kind = this->parameters_[i].kind;
        const auto j_kind = this->parameters_[j].kind;
        if(i_kind == bead_kind::Protein && j_kind == bead_kind::DNA)
        {
            return true; // Protein-DNA
        }
        else if(i_kind == bead_kind::DNA && j_kind == bead_kind::Protein)
        {
            return true; // DNA-Protein
        }
        return false; // other pair
    }

    // for testing
    exclusion_list_type const& exclusion_list() const noexcept
    {
        return exclusion_list_;
    }

    // ------------------------------------------------------------------------
    // used by Observer.
    static const char* name() noexcept {return "PDNS";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    // access to the parameters...
    std::vector<parameter_type>&       parameters()       noexcept {return parameters_;}
    std::vector<parameter_type> const& parameters() const noexcept {return parameters_;}

    std::vector<std::size_t> const& participants() const noexcept {return participants_;}
    std::vector<std::size_t> const& proteins()     const noexcept {return proteins_;}

  private:

    real_type sigma_, delta_, delta2_, one_over_2delta_, cutoff_ratio_;
    container_type           parameters_;     // indices
    std::vector<std::size_t> participants_;   // PRO + DNA beads
    std::vector<std::size_t> proteins_;       // PRO only; to reduce loop size
    exclusion_list_type      exclusion_list_;
};

} // mjolnir
#endif// MJOLNIR_FORCEFIELD_PDNS_PROTEIN_DNA_NON_SPECIFIC_INTERACTION_HPP
