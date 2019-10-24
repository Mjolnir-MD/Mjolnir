#ifndef MJOLNIR_FORCEFIELD_PWMCOS_PWMCOS_POTENTIAL_HPP
#define MJOLNIR_FORCEFIELD_PWMCOS_PWMCOS_POTENTIAL_HPP
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ExclusionList.hpp>
#include <cstdint>

namespace mjolnir
{

// Protein-DNA PWM-based interaction that represents sequence specificity of a
// protein.
// This is an implementation of the potential developed in the following paper.
// - C.Tan, and S.Takada (2018) JCTC
//
// ```toml
// [[forcefields.global]]
// interaction = "PWMcos"
// potential   = "PWMcos"
// spatial_partition.type   = "CellList"
// spatial_partition.margin = 1.0
// sigma        = 1.0     # distance sensitivity
// phi          = 0.17453 # angle sensitivity
// energy_unit  = 0.593   # overall energy coefficient (normally kBT)
// energy_shift = 0.0     # overall energy shift
// cutoff       = 5.0     # 5sigma
// parameters  = [
// {index =    0, kind = "DNA", S = 1,         B5 = 3, base = "A"},
// {index =    3, kind = "DNA", S = 4, B3 = 0, B5 = 6, base = "T"},
// # ...
// {index = 1000, kind = "Protein", CaN =  999, CaC = 1001, gamma = 1.2, epsilon = -0.4, r0 = 5.0, theta1_0 = 1.57, theta2_0 = 1.57, theta3_0 = 3.14, A = 0.5, C = 0.1, G = 0.4, T = 0.2},
// {index = 1023, kind = "Protein", CaN = 1022, CaC = 1024, gamma = 1.2, epsilon = -0.4, r0 = 6.0, theta1_0 = 1.57, theta2_0 = 1.57, theta3_0 = 3.14, A = 0.2, C = 0.4, G = 0.1, T = 0.4},
// # ...
// ]
// ```

namespace parameter_PWMcos
{
enum class base_kind : std::uint8_t
{
    A = 0,
    C = 1,
    G = 2,
    T = 3,
    X = 255,
};
} // parameter_PWMcos

template<typename traitsT>
class ProteinDNANonSpecificPotential
{
  public:
    using traits_type = traitsT;
    using real_type   = typename traits_type::real_type;
    using self_type   = ProteinDNANonSpecificPotential<traits_type>;
    using base_kind   = parameter_PWMcos::base_kind;

    struct contact_parameter_type
    {
        std::uint32_t Ca, CaN, CaC;
        real_type gamma, epsilon, r0, theta1_0, theta2_0, theta3_0;
        real_type r_cut_sq;
        std::array<real_type, 4> PWM;
    };
    struct dna_parameter_type
    {
        base_kind     base;
        std::uint32_t B, S, B5, B3;

        dna_parameter_type() : base(base_kind::X),
            B(invalid()), S(invalid()), B5(invalid()), B3(invalid())
        {}
    };
    struct parameter_type
    {
        std::uint32_t dna_index;
    };
    struct pair_parameter_type
    {
        base_kind     base;
        std::uint32_t S, B5, B3;

        pair_parameter_type(): base(base_kind::X),
            S(invalid()), B5(invalid()), B3(invalid())
        {}
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
        return parameter_type{invalid()};
    }
    static constexpr real_type default_cutoff() noexcept {return 5.0;}

  public:

    ProteinDNANonSpecificPotential(const real_type sigma, const real_type phi,
        const real_type energy_unit, const real_type energy_shift,
        const real_type cutoff_ratio,
        const std::vector<contact_parameter_type>&         contacts,
        const std::vector<dna_parameter_type>&             dnas,
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : sigma_(sigma), phi_(phi), phi2_(phi * 2),
          pi_over_2phi_(math::constants<real_type>::pi() * 0.5 / phi),
          energy_unit_(energy_unit),   energy_shift_(energy_shift),
          cutoff_ratio_(cutoff_ratio), max_cutoff_length_(0),
          contacts_(contacts),         dna_params_(dnas),
          exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {
        // this->contacts is already initialized
        this->participants_.reserve(dnas.size() + contacts.size());
        this->parameters_  .reserve(dnas.size() + contacts.size());
        this->dnas_        .reserve(dnas.size());
        this->proteins_    .reserve(contacts.size());

        for(std::size_t idx=0; idx<dnas.size(); ++idx)
        {
            const auto& dna = dnas[idx];
            this->participants_.push_back(dna.B);
            this->dnas_        .push_back(dna.B);
            if(this->parameters_.size() <= dna.B)
            {
                this->parameters_.resize(dna.B+1, default_parameter());
            }
            this->parameters_.at(dna.B).dna_index = idx;
        }
        for(const auto& contact : contacts)
        {
            this->participants_.push_back(contact.Ca);
            this->proteins_    .push_back(contact.Ca);
            if(this->parameters_.size() <= contact.Ca)
            {
                this->parameters_.resize(contact.Ca+1, default_parameter());
            }
            this->parameters_.at(contact.Ca).dna_index = invalid();
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
            const auto r_cut = para.r0 + cutoff_ratio_ * sigma_;
            para.r_cut_sq = para.r_cut * para.r_cut;
            this->max_cutoff_length_ = std::max(max_cutoff_length_, r_cut);
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
        if(abs_dtheta < this->phi_)
        {
            return std::make_pair(real_type(1), real_type(0));
        }
        else if(abs_dtheta < this->phi2_)
        {
            const real_type c = std::cos(dtheta * pi_over_2phi_);
            const real_type s = std::sin(dtheta * pi_over_2phi_);
            return std::make_pair(1 - c * c, 2 * s * c * pi_over_2phi_);
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
        if(abs_dtheta < this->phi_)
        {
            return 1;
        }
        else if(abs_dtheta < this->phi2_)
        {
            const real_type term = std::cos(dtheta * pi_over_2phi_);
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
        // should be {protein, dna}
        assert(this->parameters_.at(i).dna_index == invalid());
        assert(this->parameters_.at(j).dna_index != invalid());
        const auto& d = this->dna_params_[this->parameters_.at(j).dna_index];
        return pair_parameter_type{d.base, d.S, d.B5, d.B3};
    }

    // {PRO-Ca} U {DNA-B}
    std::vector<std::size_t> const& participants() const noexcept
    {
        return participants_;
    }
    // {Pro-Ca}
    std::vector<std::size_t> const& leading_participants() const noexcept
    {
        return this->proteins_;
    }
    // {DNA-B}
    std::vector<std::size_t> const&
    possible_partners_of(const std::size_t /*participant_idx*/,
                         const std::size_t /*particle_idx*/) const noexcept
    {
        return this->dnas_;
    }

    // to check bases has base-pairing interaction.
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept
    {
        const bool i_is_dna = (parameters_[i].dna_index != invalid());
        const bool j_is_dna = (parameters_[j].dna_index != invalid());
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
    static const char* name() noexcept {return "PWMcos";}

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

    real_type sigma_, phi_, phi2_, pi_over_2phi_;
    real_type energy_unit_, energy_shift_;
    real_type cutoff_ratio_, max_cutoff_length_;
    std::vector<std::size_t>            participants_; // all participants
    std::vector<std::size_t>            proteins_;     // indices of protein
    std::vector<std::size_t>            dnas_;         // indices of DNA
    std::vector<contact_parameter_type> contacts_;     // contact parameters
    std::vector<dna_parameter_type>     dna_params_;   // dna parameters
    std::vector<parameter_type>         parameters_;   // is_dna
    exclusion_list_type                 exclusion_list_;
};
} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize major use-cases
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

namespace mjolnir
{
extern template class PWMcosPotential<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class PWMcosPotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class PWMcosPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class PWMcosPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif// MJOLNIR_FORCEFIELD_PWMCOS_PWMCOS_POTENTIAL_HPP
