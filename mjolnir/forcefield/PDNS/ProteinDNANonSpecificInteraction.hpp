#ifndef MJOLNIR_FORCEFIELD_PDNS_PROTEIN_DNA_NON_SPECIFIC_INTERACTION_HPP
#define MJOLNIR_FORCEFIELD_PDNS_PROTEIN_DNA_NON_SPECIFIC_INTERACTION_HPP
#include <mjolnir/forcefield/PDNS/ProteinDNANonSpecificPotential.hpp>
#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/range.hpp>
#include <memory>

namespace mjolnir
{

// Protein-DNA Non-Specific interaction that represents hydrogen bond between
// side chain atoms in proteins and the phosphate residues in DNA.
// This is an implementation of the potential developed in the following paper.
// - T.Niina, G.B.Brandani, C.Tan, and S.Takada (2017) PLoS. Comp. Biol.
//
// XXX: It inherits LocalInteraction, but it is a bit different from other
// XXX: local interactions because a protein particle can have contacts with
// XXX: multiple DNA (phosphate) particles. But it also differs from global
// XXX: interactions because protein particle can have multiple interacting
// XXX: angles. The nature of this interaction is closer to a local interaction
// XXX: than a global, so it is classified as a local interaction.
//
// U(r, theta, phi) = k f(r) g(theta) g(phi)
//
// f(r)   = exp(-(r-r0)^2 / 2sigma^2)
// g(phi) = 1                                ...         |phi - phi0| <  delta
//          1 - cos^2(pi(phi-phi0) / 2delta) ... delta < |phi - phi0| < 2delta
//          0                                ... otherwise
//
template<typename traitsT>
class ProteinDNANonSpecificInteraction final : public LocalInteractionBase<traitsT>
{
  public:

    using traits_type          = traitsT;
    using base_type            = LocalInteractionBase<traits_type>;
    using real_type            = typename base_type::real_type;
    using coordinate_type      = typename base_type::coordinate_type;
    using system_type          = typename base_type::system_type;
    using topology_type        = typename base_type::topology_type;
    using connection_kind_type = typename base_type::connection_kind_type;
    using boundary_type        = typename base_type::boundary_type;
    using potential_type       = ProteinDNANonSpecificPotential<real_type>;

    // neighbor list stuff
    using neighbor_type = std::pair<std::size_t, std::size_t>; // {DNA, S3}
    using neighbor_container_type = std::vector<neighbor_type>;
    using range_type = range<typename neighbor_container_type::const_iterator>;

  public:
    ProteinDNANonSpecificInteraction(
            potential_type&& pot, const real_type margin)
        : margin_(margin), current_margin_(0), potential_(std::move(pot))
    {}
    ~ProteinDNANonSpecificInteraction() {}

    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is PDNS");

        this->potential_.initialize(sys);
        this->make_list(sys);
        return;
    }

    void update(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is PDNS");

        this->potential_.update(sys);
        this->make_list(sys);
        return;
    }

    void update_margin(const real_type dmargin, const system_type& sys) override
    {
        this->current_margin_ -= dmargin;
        if(this->current_margin_ < 0)
        {
            this->make_list(sys);
        }
        return ;
    }

    void write_topology(topology_type&) const
    {
        // no fixed topology here because contactee changes every step...
        return;
    }

    void      calc_force (system_type& sys)       const noexcept override;
    real_type calc_energy(const system_type& sys) const noexcept override;

    std::string name() const override {return "PDNSInteraction";}

    potential_type const& potential() const noexcept {return potential_;}
    potential_type&       potential()       noexcept {return potential_;}

  private:

    // TODO: it uses brute force search because there are not so many particles
    //       that have contacts. In future, Cell-list or some kind of spatial
    //       partition method would be needed ?
    void make_list(const system_type& sys)
    {
        // absolute length of margin (this->margin_ is a relative length).
        const real_type abs_margin =
            this->potential_.max_cutoff_length() * this->margin_;

        this->neighbors_.clear();
        this->neighbors_.reserve(this->potential_.contacts().size()*2);

        this->ranges_.clear();
        this->ranges_.resize(this->potential_.contacts().size()+1);
        this->ranges_.front() = 0u;
        for(std::size_t i=0; i < this->potential_.contacts().size(); ++i)
        {
            const auto& para = potential_.contacts()[i];
            const auto  P    = para.P;
            const auto& rP   = sys.position(P);
            const auto  rc   = para.r_cut + abs_margin;
            for(const auto DS : potential_.dnas())
            {
                const auto   D  = DS.first;
                const auto& rD  = sys.position(D);
                const auto dpos = sys.adjust_direction(rP - rD);
                const auto len2 = math::length_sq(dpos);

                if(len2 < rc * rc)
                {
                    this->neighbors_.push_back(DS);
                }
            }
            // assign a range for i-th contact
            this->ranges_[i+1] = this->neighbors_.size();
        }
        this->current_margin_ = abs_margin;
        return;
    }

    range_type partners_of(const std::size_t i) const noexcept
    {
        return range_type{
            this->neighbors_.begin() + this->ranges_[i],
            this->neighbors_.begin() + this->ranges_[i+1]
        };
    }

  private:

    real_type      margin_, current_margin_;
    potential_type potential_;

    neighbor_container_type  neighbors_; // dna idxs for each contacts
    std::vector<std::size_t> ranges_;    // ranges of `neighbors_`
};

template<typename traitsT>
void ProteinDNANonSpecificInteraction<traitsT>::calc_force(
        system_type& sys) const noexcept
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();

    constexpr auto tolerance = math::abs_tolerance<real_type>();
    // XXX Note: P is ambiguous because both Protein and Phosphate has `P`.
    // But this interaction is named as P-D ns, so here it uses `P` for protein
    // and `D` for DNA.

    for(std::size_t i=0; i < this->potential_.contacts().size(); ++i)
    {
        const auto& para = potential_.contacts()[i];

        const auto  P  = para.P;
        const auto& rP = sys.position(P);
        for(const auto& ptnr : this->partners_of(i))
        {
            const auto  D  = ptnr.first;  // DNA phosphate
            const auto  S3 = ptnr.second; // 3' Sugar
            const auto& rD = sys.position(D);

            MJOLNIR_LOG_DEBUG("protein = ", P, ", DNA = ", D, ", r0 = ", para.r0);

            //  PC          S5'    |
            //    o         o--o B | theta is an angle formed by the vector
            //    A\ P   D /       | from the neighboring amino acid residue at
            // vP | o --> o        | the N-term side to the amino acid residue
            //    |/     `-\       | at the C-term side (vP in the figure) and
            //    o    phi  o--o B | the vector from Protein to DNA (phosphate).
            //  PN           S3'   |

            // ----------------------------------------------------------------
            // calculates the distance part

            const auto rPD    = sys.adjust_direction(rD - rP); // PRO -> DNA
            const auto lPD_sq = math::length_sq(rPD);
            if(para.r_cut_sq < lPD_sq)
            {
                continue;
            }
            const auto rlPD   = math::rsqrt(lPD_sq);
            const auto  lPD   = lPD_sq * rlPD;
            const auto f_df   = potential_.f_df(para.r0, lPD);

            MJOLNIR_LOG_DEBUG("f = ", f_df.first, ", df = ", f_df.second);

            // ----------------------------------------------------------------
            // calculates the angle part (theta)

            const auto& rPC   = sys.position(para.PC);
            const auto& rPN   = sys.position(para.PN);
            const auto rPNC   = sys.adjust_direction(rPC - rPN); // PN -> PC
            const auto rlPNC  = math::rlength(rPNC);
            const auto dotPNC = math::dot_product(rPNC, rPD);
            const auto cosPNC = dotPNC * rlPD * rlPNC;
            const auto theta  = std::acos(math::clamp<real_type>(cosPNC,-1,1));

            const auto g_dg_theta = potential_.g_dg(para.theta0, theta);

            MJOLNIR_LOG_DEBUG("g(theta) = ", g_dg_theta.first,
                             ", dg(theta) = ", g_dg_theta.second);

            // ----------------------------------------------------------------
            // calculates the angle part (phi)

            const auto& rS3   = sys.position(S3);
            const auto rS3D   = sys.adjust_direction(rD - rS3); // S3' -> D
            const auto rlS3D  = math::rlength(rS3D);
            const auto dotS3D = math::dot_product(rPD, rS3D);
            const auto cosS3D = dotS3D * rlS3D * rlPD;
            const auto phi    = std::acos(math::clamp<real_type>(cosS3D,-1,1));

            const auto g_dg_phi = potential_.g_dg(para.phi0, phi);

            MJOLNIR_LOG_DEBUG("g(phi) = ", g_dg_phi.first,
                             ", dg(phi) = ", g_dg_phi.second);

            // ----------------------------------------------------------------
            // calculate force
            //
            //   d/dr [kf(r)  g(theta)  g(phi)]
            // =    k [df(r)  g(theta)  g(phi) dr/dr     +
            //          f(r) dg(theta)  g(phi) dtheta/dr +
            //          f(r)  g(theta) dg(phi) dphi/dr   ]
            const auto k = para.k;

            // df(r) g(theta) g(phi)
            if(g_dg_theta.first  != 0 && g_dg_phi.first  != 0)
            {
                MJOLNIR_LOG_DEBUG("calculating distance force");

                const auto coef = rlPD * k *
                    f_df.second * g_dg_theta.first * g_dg_phi.first;
                const auto F = -coef * rPD;

                sys.force(P) -= F;
                sys.force(D) += F;
            }

            // f(r) dg(theta) g(phi)
            if(g_dg_theta.second != 0 && g_dg_phi.first  != 0)
            {
                MJOLNIR_LOG_DEBUG("calculating theta force");

                const auto deriv =
                    k * f_df.first * g_dg_theta.second * g_dg_phi.first;

                const auto sin_theta = std::sin(theta);
                const auto coef_sin  = deriv / std::max(sin_theta, tolerance);

                const auto rPD_reg  = rlPD  * rPD;
                const auto rPNC_reg = rlPNC * rPNC;

                const auto F_P  = -coef_sin * rlPD  * (rPNC_reg - cosPNC * rPD_reg );
                const auto F_PN = -coef_sin * rlPNC * (rPD_reg  - cosPNC * rPNC_reg);

                sys.force(D)       -= F_P;
                sys.force(P)       += F_P;
                sys.force(para.PN) += F_PN;
                sys.force(para.PC) -= F_PN;
            }

            // f(r) dg(theta) g(phi)
            if(g_dg_theta.first  != 0 && g_dg_phi.second != 0)
            {
                MJOLNIR_LOG_DEBUG("calculating phi force");

                const auto deriv =
                    k * f_df.first * g_dg_theta.first * g_dg_phi.second;

                const auto sin_phi  = std::sin(phi);
                const auto coef_sin = deriv / std::max(sin_phi, tolerance);

                const auto rPD_reg  = rlPD  * rPD;
                const auto rS3D_reg = rlS3D * rS3D;

                const auto F_P = -coef_sin * rlPD  * (rS3D_reg - cosS3D * rPD_reg);
                const auto F_S = -coef_sin * rlS3D * (rPD_reg  - cosS3D * rS3D_reg);

                sys.force(P)  += F_P;
                sys.force(D)  -= F_P + F_S;
                sys.force(S3) += F_S;
            }
        }
    }
    return;
}

template<typename traitsT>
typename ProteinDNANonSpecificInteraction<traitsT>::real_type
ProteinDNANonSpecificInteraction<traitsT>::calc_energy(
        const system_type& sys) const noexcept
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    // XXX Note: P is ambiguous because both Protein and Phosphate has `P`.
    // But this interaction is named as P-D ns, so here it uses `P` for protein
    // and `D` for DNA.

    real_type E = 0.0;
    for(std::size_t i=0; i < this->potential_.contacts().size(); ++i)
    {
        const auto& para = potential_.contacts()[i];

        const auto  P  = para.P;
        const auto& rP = sys.position(P);
        for(const auto& ptnr : this->partners_of(i))
        {
            const auto  D  = ptnr.first;
            const auto  S3 = ptnr.second;
            const auto& rD = sys.position(D);

            //  PC          S5'    |
            //    o         o--o B | theta is an angle formed by the vector
            //    A\ P   D /       | from the neighboring amino acid residue at
            // vP | o --> o        | the N-term side to the amino acid residue
            //    |/     `-\       | at the C-term side (vP in the figure) and
            //    o    phi  o--o B | the vector from Protein to DNA (phosphate).
            //  PN           S3'   |

            // ----------------------------------------------------------------
            // calculates the distance part

            const auto rPD    = sys.adjust_direction(rD - rP); // PRO -> DNA
            const auto lPD_sq = math::length_sq(rPD);
            if(para.r_cut_sq < lPD_sq)
            {
                continue;
            }
            const auto rlPD   = math::rsqrt(lPD_sq);
            const auto  lPD   = lPD_sq * rlPD;
            const auto f      = potential_.f(para.r0, lPD);

            // ----------------------------------------------------------------
            // calculates the angle part (theta)

            const auto& rPC    = sys.position(para.PC);
            const auto& rPN    = sys.position(para.PN);
            const auto rPNC    = sys.adjust_direction(rPC - rPN); // PN -> PC
            const auto rlPNC   = math::rlength(rPNC);
            const auto dotPNC  = math::dot_product(rPNC, rPD);
            const auto cosPNC  = dotPNC * rlPD * rlPNC;
            const auto theta   = std::acos(math::clamp<real_type>(cosPNC,-1,1));
            const auto g_theta = potential_.g(para.theta0, theta);

            if(g_theta == real_type(0)) {continue;}

            // ----------------------------------------------------------------
            // calculates the angle part (phi)

            const auto& rS3   = sys.position(S3);
            const auto rS3D   = sys.adjust_direction(rD - rS3); // S3' -> D
            const auto rlS3D  = math::rlength(rS3D);
            const auto dotS3D = math::dot_product(rPD, rS3D);
            const auto cosS3D = dotS3D * rlS3D * rlPD;
            const auto phi    = std::acos(math::clamp<real_type>(cosS3D,-1,1));
            const auto g_phi  = potential_.g(para.phi0, phi);

            if(g_phi == real_type(0)) {continue;}

            // ----------------------------------------------------------------
            // calculate energy

            const auto k = para.k;

            MJOLNIR_LOG_INFO("protein = ", P, ", DNA = ", D, ", r0 = ", para.r0);

            E += k * f * g_theta * g_phi;
        }
    }
    return E;
}

} // mjolnir
#endif// MJOLNIR_FORCEFIELD_PDNS_PROTEIN_DNA_NON_SPECIFIC_INTERACTION_HPP
