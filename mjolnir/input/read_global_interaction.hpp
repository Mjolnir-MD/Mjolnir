#ifndef MJOLNIR_INPUT_READ_GLOBAL_INTERACTION_HPP
#define MJOLNIR_INPUT_READ_GLOBAL_INTERACTION_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/forcefield/global/GlobalPairInteraction.hpp>
#include <mjolnir/forcefield/global/GlobalPairLennardJonesInteraction.hpp>
#include <mjolnir/forcefield/global/GlobalPairUniformLennardJonesInteraction.hpp>
#include <mjolnir/forcefield/global/GlobalPairExcludedVolumeInteraction.hpp>
#include <mjolnir/forcefield/3SPN2/ThreeSPN2BaseBaseInteraction.hpp>
#include <mjolnir/forcefield/PDNS/ProteinDNANonSpecificInteraction.hpp>
#include <mjolnir/forcefield/PWMcos/PWMcosInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_global_potential.hpp>
#include <mjolnir/input/read_spatial_partition.hpp>
#include <memory>

namespace mjolnir
{

// ----------------------------------------------------------------------------
// global interaction
// ----------------------------------------------------------------------------

template<typename traitsT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_pair_interaction(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto potential = toml::find<std::string>(global, "potential");

    if(potential == "ExcludedVolume")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Excluded Volume.");
        using potential_t   = ExcludedVolumePotential<traitsT>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        return make_unique<interaction_t>(
            read_excluded_volume_potential<traitsT>(global),
            read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "InversePower")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Inverse Power.");
        using potential_t   = InversePowerPotential<traitsT>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        return make_unique<interaction_t>(
            read_inverse_power_potential<traitsT>(global),
            read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "HardCoreExcludedVolume")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Hard Core Excluded Volume.");
        using potential_t   = HardCoreExcludedVolumePotential<traitsT>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        return make_unique<interaction_t>(
            read_hard_core_excluded_volume_potential<traitsT>(global),
            read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "DebyeHuckel")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Debye-Huckel.");
        using potential_t   = DebyeHuckelPotential<traitsT>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        return make_unique<interaction_t>(
            read_debye_huckel_potential<traitsT>(global),
            read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "LennardJones")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Lennard-Jones.");
        using potential_t   = LennardJonesPotential<traitsT>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        return make_unique<interaction_t>(
            read_lennard_jones_potential<traitsT>(global),
            read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "UniformLennardJones")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Uniform Lennard-Jones.");
        using potential_t   = UniformLennardJonesPotential<traitsT>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        return make_unique<interaction_t>(
            read_uniform_lennard_jones_potential<traitsT>(global),
            read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "3SPN2ExcludedVolume")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is 3SPN2ExcludedVolume.");
        using potential_t   = ThreeSPN2ExcludedVolumePotential<traitsT>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        return make_unique<interaction_t>(
            read_3spn2_excluded_volume_potential<traitsT>(global),
            read_spatial_partition<traitsT, potential_t>(global));
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_global_pair_interaction: invalid potential",
            toml::find<toml::value>(global, "potential"), "here", {
            "expected value is one of the following.",
            "- \"ExcludedVolume\"       : repulsive r^12 potential",
            "- \"InversePower\"         : repulsive r^n potential",
            "- \"DebyeHuckel\"          : Debye-Huckel type electrostatic potential",
            "- \"LennardJones\"         : famous r^12 - r^6 potential",
            "- \"UniformLennardJones\"  : famous r^12 - r^6 potential with uniform parameters",
            "- \"3SPN2ExcludedVolume\"  : excluded volume for 3SPN2 DNA model"
            }));
    }
}

// ----------------------------------------------------------------------------
// 3SPN2 Base-Base Interaction

template<typename traitsT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_3spn2_base_base_interaction(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type           = typename traitsT::real_type;
    using base_kind           = parameter_3SPN2::base_kind;
    using potential_type      = ThreeSPN2BaseBaseInteractionPotential<traitsT>;
    using parameter_type      = typename potential_type::parameter_type;

    // [[forcefields.global]]
    // interaction = "3SPN2BaseBase"
    // potential   = "3SPN2"
    // spatial_partition = {type = "CellList", margin = 1.0}
    // parameters = [
    // {strand = 0, nucleotide =  0,          S =   0, B =   1, Base = "A"},
    // {strand = 0, nucleotide =  1, P =   2, S =   3, B =   4, Base = "T"},
    // # ...
    // ]

    // ------------------------------------------------------------------------
    // read parameters

    const auto& env    = global.as_table().count("env") == 1 ?
                         global.as_table().at("env") : toml::value{};

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    using nucleotide_index_type = parameter_3SPN2::NucleotideIndex;
    std::vector<nucleotide_index_type> nuc_idxs;

    nuc_idxs.reserve(ps.size());
    for(const auto& item : ps)
    {
        nucleotide_index_type nuc_idx;
        const auto ofs = find_parameter_or<std::int64_t>(item, env, "offset", 0);

        // at the edge of the DNA, Phosphate may not exist.
        if(item.as_table().count("P") != 0)
        {
            nuc_idx.P = find_parameter<std::size_t>(item, env, "P") + ofs;
        }
        nuc_idx.S          = find_parameter<std::size_t>(item, env, "S") + ofs;
        nuc_idx.B          = find_parameter<std::size_t>(item, env, "B") + ofs;
        nuc_idx.strand     = find_parameter<std::size_t>(item, env, "strand");
        nuc_idx.nucleotide = find_parameter<std::size_t>(item, env, "nucleotide");

        const auto bk      = toml::find<std::string>(item, "Base");
        if     (bk == "A") {nuc_idx.base = base_kind::A;}
        else if(bk == "T") {nuc_idx.base = base_kind::T;}
        else if(bk == "G") {nuc_idx.base = base_kind::G;}
        else if(bk == "C") {nuc_idx.base = base_kind::C;}
        else
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_local_3spn2_base_base_interaction: "
                "invalid Base", item, "here", {
                "expected value is one of the \"A\", \"T\", \"C\", \"G\"."
                }));
        }

        MJOLNIR_LOG_INFO("ThreeSPN2BaseBaseInteraction: nucleotide = ", nuc_idx);
        nuc_idxs.push_back(nuc_idx);
    }

    std::sort(nuc_idxs.begin(), nuc_idxs.end(),
        [](const nucleotide_index_type& lhs, const nucleotide_index_type& rhs) {
            return std::make_pair(lhs.strand, lhs.nucleotide) <
                   std::make_pair(rhs.strand, rhs.nucleotide);
        });

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(nuc_idxs.size());
    for(std::size_t i=0; i<nuc_idxs.size(); ++i)
    {
        const auto& base = nuc_idxs.at(i);
        const auto  B = base.B;

        parameter_type p;
        p.nucleotide_index = base.nucleotide;
        p.base   = base.base;
        p.S_idx  = base.S;
        p.B5_idx = potential_type::invalid();
        p.B3_idx = potential_type::invalid();

        if(i != 0 && nuc_idxs.at(i-1).strand == base.strand)
        {
            p.B5_idx = nuc_idxs.at(i-1).B;
        }
        if(i+1 < nuc_idxs.size() && nuc_idxs.at(i+1).strand == base.strand)
        {
            p.B3_idx = nuc_idxs.at(i+1).B;
        }
        MJOLNIR_LOG_INFO("Base idx = ", B, ", base = ", p.base, ", Sugar idx = ",
            p.S_idx, ", 5' adjacent = ", p.B5_idx, ", 3' adjacent", p.B3_idx);

        params.emplace_back(B, p);
    }

    const auto pot = toml::find<std::string>(global, "potential");
    if(pot == "3SPN2")
    {
        ThreeSPN2BaseBaseGlobalPotentialParameter<real_type> para_3SPN2;
        potential_type potential(para_3SPN2, std::move(params),
            read_ignore_particles_within(global),
            read_ignored_molecule(global), read_ignored_group(global));

        return make_unique<ThreeSPN2BaseBaseInteraction<traitsT>>(
                std::move(potential),
                read_spatial_partition<traitsT, potential_type>(global));
    }
    else if(pot == "3SPN2C")
    {
        ThreeSPN2CBaseBaseGlobalPotentialParameter<real_type> para_3SPN2C;
        potential_type potential(para_3SPN2C, std::move(params),
            read_ignore_particles_within(global),
            read_ignored_molecule(global), read_ignored_group(global));

        return make_unique<ThreeSPN2BaseBaseInteraction<traitsT>>(
                std::move(potential),
                read_spatial_partition<traitsT, potential_type>(global));
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_local_3spn2_base_stacking_interaction: "
            "invalid potential", toml::find(global, "potential"), "here", {
            "expected value is one of the following.",
            "- \"3SPN2\" : The general 3SPN2 parameter set.",
            "- \"3SPN2C\": The parameter set optimized to reproduce curveture of dsDNA."
            }));
    }
}

// ----------------------------------------------------------------------------
// PDNS Interaction

template<typename traitsT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_pdns_interaction(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using potential_type         = ProteinDNANonSpecificPotential<traitsT>;
    using real_type              = typename potential_type::real_type;
    using contact_parameter_type = typename potential_type::contact_parameter_type;
    using dna_index_type         = typename potential_type::dna_index_type;

    // ```toml
    // [[forcefields.global]]
    // interaction = "PDNS"
    // potential   = "PDNS"
    // spatial_partition.type = "VerletList"
    // spatial_partition.margin = 0.4
    // sigma  = 1.0
    // delta  = 0.17453
    // cutoff = 5.0 # relative to sigma
    // parameters  = [
    // {index =    2, kind = "DNA", S3 = 1},
    // {index =    5, kind = "DNA", S3 = 4},
    // # ...
    // {index = 1000, kind = "Protein", PN =  999, PC = 1001, k = 1.2, r0 = 5.0, theta0 = 100.0, phi0 = 130.0},
    // {index = 1023, kind = "Protein", PN = 1022, PC = 1024, k = 1.2, r0 = 6.0, theta0 = 110.0, phi0 = 120.0},
    // # ...
    // ]
    // ```

    // ------------------------------------------------------------------------
    // read parameters

    const real_type dlt = toml::find<real_type>(global, "delta");
    const real_type sgm = toml::find<real_type>(global, "sigma");

    const real_type cutoff = toml::find_or<real_type>(global, "cutoff",
            potential_type::default_cutoff());

    const auto& env = global.as_table().count("env") == 1 ?
                      global.as_table().at("env") : toml::value{};

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<contact_parameter_type> contacts;
    std::vector<dna_index_type>         dnas;
    contacts.reserve(ps.size());
    dnas    .reserve(ps.size());

    for(const auto& item : ps)
    {
        const auto idx  = toml::find<std::size_t>(item, "index");
        const auto ofs  = toml::find_or<std::int64_t>(item, "offset", 0);
        const auto kind = toml::find<std::string>(item, "kind");
        if(kind == "Protein")
        {
            contact_parameter_type para;
            para.P      = idx + ofs;
            para.PN     = find_parameter<std::uint32_t>(item, env, "PN") + ofs;
            para.PC     = find_parameter<std::uint32_t>(item, env, "PC") + ofs;
            para.k      = find_parameter<real_type>(item, env, "k");
            para.r0     = find_parameter<real_type>(item, env, "r0");
            para.theta0 = find_parameter<real_type>(item, env, "theta0");
            para.phi0   = find_parameter<real_type>(item, env, "phi0");
            contacts.push_back(para);

            MJOLNIR_LOG_INFO("Protein: idx = ", para.P, ", PN = ", para.PN,
                ", PC = ", para.PC, ", k = ", para.k, ", r0 = ", para.r0,
                ", theta0 = ", para.theta0, ", phi0 = ", para.phi0);
        }
        else if (kind == "DNA")
        {
            dna_index_type di;
            di.D  = idx + ofs;
            di.S3 = find_parameter<std::uint32_t>(item, env, "S3") + ofs;
            dnas.push_back(di);
            MJOLNIR_LOG_INFO("DNA: idx = ", di.D, ", S3 = ", di.S3);
        }
        else
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_pdns_interaction: unknown kind ",
                toml::find(item, "kind"), "here", {
                "expected value is one of the following.",
                "- \"Protein\": Protein bead has PN, PC and native parameters",
                "- \"DNA\"    : DNA bead has index of the corresponding Sugar."
                }));
        }
    }
    return make_unique<ProteinDNANonSpecificInteraction<traitsT>>(
        potential_type(sgm, dlt, cutoff, std::move(contacts), std::move(dnas),
            read_ignore_particles_within(global), read_ignored_molecule(global),
            read_ignored_group(global)),
        read_spatial_partition<traitsT, potential_type>(global));
}

// ----------------------------------------------------------------------------
// PWMcos Interaction

template<typename traitsT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_pwmcos_interaction(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using potential_type         = PWMcosPotential<traitsT>;
    using real_type              = typename potential_type::real_type;
    using contact_parameter_type = typename potential_type::contact_parameter_type;
    using dna_parameter_type     = typename potential_type::dna_parameter_type;
    using base_kind              = typename potential_type::base_kind;

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
    // {index =    0, S = 1,         B5 = 3, base = "A"},
    // {index =    3, S = 4, B3 = 0, B5 = 6, base = "T"},
    // # ...
    // {index = 1000, PN =  999, PC = 1001, gamma = 1.2, epsilon = -0.4, r0 = 5.0, theta1_0 = 1.57, theta2_0 = 1.57, theta3_0 = 3.14, A = 0.5, C = 0.1, G = 0.4, T = 0.2},
    // {index = 1023, PN = 1022, PC = 1024, gamma = 1.2, epsilon = -0.4, r0 = 6.0, theta1_0 = 1.57, theta2_0 = 1.57, theta3_0 = 3.14, A = 0.2, C = 0.4, G = 0.1, T = 0.4},
    // # ...
    // ]
    // ```

    // ------------------------------------------------------------------------
    // read parameters

    const real_type phi    = toml::find<real_type>(global, "phi");
    const real_type sgm    = toml::find<real_type>(global, "sigma");
    const real_type Eunit  = toml::find<real_type>(global, "energy_unit");
    const real_type Eshift = toml::find<real_type>(global, "energy_shift");

    const real_type cutoff = toml::find_or<real_type>(global, "cutoff",
            potential_type::default_cutoff());

    const auto& env = global.as_table().count("env") == 1 ?
                      global.as_table().at("env") : toml::value{};

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<contact_parameter_type> contacts;
    std::vector<dna_parameter_type>     dnas;
    contacts.reserve(ps.size());
    dnas    .reserve(ps.size());

    // base index
    const auto A = static_cast<std::size_t>(base_kind::A);
    const auto C = static_cast<std::size_t>(base_kind::C);
    const auto G = static_cast<std::size_t>(base_kind::G);
    const auto T = static_cast<std::size_t>(base_kind::T);

    for(const auto& item : ps)
    {
        const auto idx  = toml::find<std::size_t>(item, "index");
        const auto ofs  = toml::find_or<std::int64_t>(item, "offset", 0);
        const auto kind = toml::find<std::string>(item, "kind");
        if(kind == "Protein")
        {
            contact_parameter_type para;
            para.Ca       = idx + ofs;
            para.CaN      = find_parameter<std::uint32_t>(item, env, "CaN") + ofs;
            para.CaC      = find_parameter<std::uint32_t>(item, env, "CaC") + ofs;
            para.gamma    = find_parameter<real_type>(item, env, "gamma");
            para.epsilon  = find_parameter<real_type>(item, env, "epsilon");
            para.r0       = find_parameter<real_type>(item, env, "r0");
            para.theta1_0 = find_parameter<real_type>(item, env, "theta1_0");
            para.theta2_0 = find_parameter<real_type>(item, env, "theta2_0");
            para.theta3_0 = find_parameter<real_type>(item, env, "theta3_0");
            para.PWM[A]   = find_parameter<real_type>(item, env, "A");
            para.PWM[C]   = find_parameter<real_type>(item, env, "C");
            para.PWM[G]   = find_parameter<real_type>(item, env, "G");
            para.PWM[T]   = find_parameter<real_type>(item, env, "T");

            contacts.push_back(para);

            MJOLNIR_LOG_INFO("Protein: idx = ", idx + ofs, ", CaN = ", para.CaN,
                ", CaC = ", para.CaC, ", gamma = ", para.gamma, ", epsilon = ", para.epsilon,
                ", r0 = ", para.r0, ", theta1_0 = ", para.theta1_0,
                ", tehta2_0 = ", para.theta2_0, ", theta3_0 = ", para.theta3_0,
                ", PWM_A = ", para.PWM[A], ", PWM_C = ", para.PWM[C],
                ", PWM_G = ", para.PWM[G], ", PWM_T = ", para.PWM[T]);
        }
        else if (kind == "DNA")
        {
            dna_parameter_type dp;
            dp.B  = idx + ofs;
            dp.S  = find_parameter<std::uint32_t>(item, env, "S") + ofs;

            // it seems that 3' and 5' bases are required.
            dp.B5 = find_parameter<std::uint32_t>(item, env, "B5") + ofs;
            dp.B3 = find_parameter<std::uint32_t>(item, env, "B3") + ofs;

            const auto& bk = toml::find<std::string>(item, "base");
            if     (bk == "A") {dp.base = base_kind::A;}
            else if(bk == "C") {dp.base = base_kind::C;}
            else if(bk == "G") {dp.base = base_kind::G;}
            else if(bk == "T") {dp.base = base_kind::T;}
            else
            {
                throw_exception<std::runtime_error>(toml::format_error("[error] "
                    "mjolnir::read_pdns_interaction: unknown base; "
                    "none of \"A\", \"T\", \"C\", \"G\"",
                    toml::find(item, "base"), "here"));
            }
            dnas.push_back(dp);
            MJOLNIR_LOG_INFO("DNA: idx = ", dp.B, ", base = ", bk,
                             ", S = ", dp.S, ", B3 = ", dp.B3, ", B5 = ", dp.B5);
        }
        else
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_pdns_interaction: unknown kind ",
                toml::find(item, "kind"), "here", {
                "expected value is one of the following.",
                "- \"Protein\": Protein bead has PN, PC and native parameters",
                "- \"DNA\"    : DNA bead has index of the corresponding Sugar."
                }));
        }
    }
    return make_unique<PWMcosInteraction<traitsT>>(
        potential_type(sgm, phi, Eunit, Eshift, cutoff,
            std::move(contacts), std::move(dnas),
            read_ignore_particles_within(global), read_ignored_molecule(global),
            read_ignored_group(global)),
        read_spatial_partition<traitsT, potential_type>(global));
}


// ----------------------------------------------------------------------------
// general read_global_interaction function
// ----------------------------------------------------------------------------

template<typename traitsT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_interaction(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    const auto interaction = toml::find<std::string>(global, "interaction");

    if(interaction == "Pair")
    {
        MJOLNIR_LOG_NOTICE("Pair interaction found.");
        return read_global_pair_interaction<traitsT>(global);
    }
    else if(interaction == "3SPN2BaseBase")
    {
        MJOLNIR_LOG_NOTICE("3SPN2BaseBaseInteraction found.");
        return read_global_3spn2_base_base_interaction<traitsT>(global);
    }
    else if(interaction == "PDNS")
    {
        MJOLNIR_LOG_NOTICE("P-D ns Interaction found.");
        return read_pdns_interaction<traitsT>(global);
    }
    else if(interaction == "PWMcos")
    {
        MJOLNIR_LOG_NOTICE("PWMcos Interaction found.");
        return read_pwmcos_interaction<traitsT>(global);
    }
    else
    {
        throw std::runtime_error(toml::format_error("[error] "
            "mjolnir::read_global_interaction: invalid interaction",
            toml::find<toml::value>(global, "interaction"), "here", {
            "expected value is one of the following.",
            "- \"Pair\"         : well-known pair interaction depends only on the distance",
            "- \"3SPN2BaseBase\": Base pair and cross stacking interaction for 3SPN2 DNA model"
            }));
    }
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_global_interaction(const toml::value& global);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_global_interaction(const toml::value& global);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_global_interaction(const toml::value& global);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_global_interaction(const toml::value& global);

extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_global_3spn2_base_base_interaction(const toml::value& global);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_global_3spn2_base_base_interaction(const toml::value& global);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_global_3spn2_base_base_interaction(const toml::value& global);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_global_3spn2_base_base_interaction(const toml::value& global);

extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_global_pair_interaction(const toml::value& global);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_global_pair_interaction(const toml::value& global);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_global_pair_interaction(const toml::value& global);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_global_pair_interaction(const toml::value& global);
#endif

} // mjolnir
#endif// MJOLNIR_READ_GLOBAL_INTERACTION_HPP
