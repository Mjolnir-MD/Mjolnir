#ifndef MJOLNIR_FORCEFIELD_PDNS_PROTEIN_DNA_NON_SPECIFIC_POTENTIAL_HPP
#define MJOLNIR_FORCEFIELD_PDNS_PROTEIN_DNA_NON_SPECIFIC_POTENTIAL_HPP
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ExclusionList.hpp>
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
//
// XXX one particle on protein can have multiple contact direction.

template<typename traitsT>
class ProteinDNANonSpecificPotential
{
  public:
    using traits_type = traitsT;
    using real_type   = typename traits_type::real_type;
    using self_type   = ProteinDNANonSpecificPotential<traits_type>;

    struct contact_parameter_type
    {
        std::uint32_t P, PN, PC;
        real_type k, r0, theta0, phi0;
        real_type r_cut, r_cut_sq;
    };
    struct dna_index_type
    {
        std::uint32_t D, S3;
    };
    struct parameter_type
    {
        std::uint32_t S3; // if a particle is a protein, set `invalid()`
    };
    struct pair_parameter_type
    {
        std::uint32_t S3;
    };
    using container_type = std::vector<parameter_type>;

    // topology stuff
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList<traits_type>;

    static constexpr std::uint32_t invalid() noexcept
    {
        return std::numeric_limits<std::uint32_t>::max();
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{invalid()};
    }
    static constexpr real_type default_cutoff() noexcept {return 5.0;}

  public:

    ProteinDNANonSpecificPotential(const real_type sigma,
        const real_type delta, const real_type cutoff_ratio,
        const std::vector<contact_parameter_type>&         contacts,
        const std::vector<dna_index_type>&                 dnas,
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : sigma_(sigma), delta_(delta), delta2_(delta * 2),
          pi_over_2delta_(math::constants<real_type>::pi() * 0.5 / delta),
          cutoff_ratio_(cutoff_ratio), max_cutoff_length_(0),
          contacts_(contacts),
          exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {
        // this->contacts is already initialized
        this->participants_.reserve(dnas.size() + contacts.size());
        this->parameters_  .reserve(dnas.size() + contacts.size());
        this->dnas_        .reserve(dnas.size());
        this->proteins_    .reserve(contacts.size());
        for(const auto& dna : dnas)
        {
            this->participants_.push_back(dna.D);
            this->dnas_        .push_back(dna.D);
            if(this->parameters_.size() <= dna.D)
            {
                this->parameters_.resize(dna.D+1, default_parameter());
            }
            this->parameters_.at(dna.D).S3 = dna.S3;
        }
        for(const auto& contact : contacts)
        {
            this->participants_.push_back(contact.P);
            this->proteins_    .push_back(contact.P);
            if(this->parameters_.size() <= contact.P)
            {
                this->parameters_.resize(contact.P+1, default_parameter());
            }
            this->parameters_.at(contact.P).S3 = invalid();
        }

        // Protein particles can have multiple contact sites, so `participants_`
        // and `proteins_` may have duplicated indices.

        std::sort(participants_.begin(), participants_.end());
        const auto uniq_participants =
            std::unique(participants_.begin(), participants_.end());
        participants_.erase(uniq_participants, participants_.end());

        std::sort(proteins_.begin(), proteins_.end());
        const auto uniq_proteins =
            std::unique(proteins_.begin(), proteins_.end());
        proteins_.erase(uniq_proteins, proteins_.end());
    }
    ~ProteinDNANonSpecificPotential() = default;

    void initialize(const System<traits_type>& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // set cutoff length

        this->max_cutoff_length_ = real_type(0);
        for(auto& para : this->contacts_)
        {
            para.r_cut    = para.r0 + cutoff_ratio_ * sigma_;
            para.r_cut_sq = para.r_cut * para.r_cut;

            this->max_cutoff_length_ = std::max(max_cutoff_length_, para.r_cut);
        }

        this->update(sys);
        return;
    }

    void update(const System<traits_type>&) noexcept
    {
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

    std::vector<contact_parameter_type> const& contacts() const noexcept
    {
        return contacts_;
    }

    // -----------------------------------------------------------------------
    // for NeighborList

    pair_parameter_type
    prepare_params(const std::size_t i, const std::size_t j) const noexcept
    {
        assert(this->parameters_.at(i).S3 == invalid());
        assert(this->parameters_.at(j).S3 != invalid());
        return pair_parameter_type{this->parameters_[j].S3};
    }

    // {PRO-Ca} U {DNA-P}
    std::vector<std::size_t> const& participants() const noexcept
    {
        return participants_;
    }
    // {Pro-Ca}
    std::vector<std::size_t> const& leading_participants() const noexcept
    {
        return this->proteins_;
    }
    // {DNA-P}
    std::vector<std::size_t> const&
    possible_partners_of(const std::size_t /*participant_idx*/,
                         const std::size_t /*particle_idx*/) const noexcept
    {
        return this->dnas_;
    }

    // to check bases has base-pairing interaction.
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept
    {
        const bool i_is_dna = (parameters_[i].S3 != invalid());
        const bool j_is_dna = (parameters_[j].S3 != invalid());
        // {protein, dna} pair has interaction, others not.
        if(!i_is_dna && j_is_dna)
        {
            return true;
        }
        assert(!i_is_dna || j_is_dna); // {dna, pro} pair should not be listed
        return false;
    }

    exclusion_list_type const& exclusion_list() const noexcept {return exclusion_list_;}

    real_type cutoff_ratio()      const noexcept {return cutoff_ratio_;}
    real_type max_cutoff_length() const noexcept {return this->max_cutoff_length_;}

    // ------------------------------------------------------------------------
    // used by Observer.
    static const char* name() noexcept {return "PDNS";}

    // ------------------------------------------------------------------------
    // for tests

    std::vector<parameter_type> const& parameters() const noexcept
    {
        return parameters_;
    }
    std::vector<std::size_t> const& dnas() const noexcept
    {
        return dnas_;
    }
    std::vector<std::size_t> const& proteins() const noexcept
    {
        return proteins_;
    }

  private:

    real_type sigma_, delta_, delta2_, pi_over_2delta_;
    real_type cutoff_ratio_, max_cutoff_length_;
    std::vector<std::size_t> participants_;          // all participants
    std::vector<std::size_t> proteins_;              // indices of protein
    std::vector<std::size_t> dnas_;                  // indices of DNA
    std::vector<contact_parameter_type> contacts_;   // contact parameters
    std::vector<parameter_type>         parameters_; // S3 index(DNA) or
                                                     // `invalid`(Protein)
    exclusion_list_type exclusion_list_;
};
} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize major use-cases
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

namespace mjolnir
{
extern template class ProteinDNANonSpecificPotential<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class ProteinDNANonSpecificPotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class ProteinDNANonSpecificPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class ProteinDNANonSpecificPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif// MJOLNIR_FORCEFIELD_PDNS_PROTEIN_DNA_NON_SPECIFIC_INTERACTION_HPP
