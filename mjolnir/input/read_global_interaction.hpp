#ifndef MJOLNIR_INPUT_READ_GLOBAL_INTERACTION_HPP
#define MJOLNIR_INPUT_READ_GLOBAL_INTERACTION_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/forcefield/global/GlobalPairInteraction.hpp>
#include <mjolnir/forcefield/global/GlobalPairLennardJonesInteraction.hpp>
#include <mjolnir/forcefield/global/GlobalPairExcludedVolumeInteraction.hpp>
#include <mjolnir/forcefield/3SPN2/ThreeSPN2BasePairInteraction.hpp>
#include <mjolnir/forcefield/3SPN2/ThreeSPN2CrossStackingInteraction.hpp>
#include <mjolnir/forcefield/PDNS/ProteinDNANonSpecificInteraction.hpp>
#include <mjolnir/forcefield/PWMcos/PWMcosInteraction.hpp>
#include <mjolnir/forcefield/stoichiometric/GlobalStoichiometricInteraction.hpp>
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
        using real_type     = typename traitsT::real_type;
        using potential_t   = ExcludedVolumePotential<real_type>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        auto pot_para = read_excluded_volume_potential<traitsT>(global);

        return make_unique<interaction_t>(
                std::move(pot_para.first), std::move(pot_para.second),
                read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "InversePower")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Inverse Power.");
        using real_type     = typename traitsT::real_type;
        using potential_t   = InversePowerPotential<real_type>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        auto pot_para = read_inverse_power_potential<traitsT>(global);

        return make_unique<interaction_t>(
                std::move(pot_para.first), std::move(pot_para.second),
                read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "HardCoreExcludedVolume")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Hard Core Excluded Volume.");
        using real_type     = typename traitsT::real_type;
        using potential_t   = HardCoreExcludedVolumePotential<real_type>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        auto pot_para = read_hard_core_excluded_volume_potential<traitsT>(global);

        return make_unique<interaction_t>(
                std::move(pot_para.first), std::move(pot_para.second),
                read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "DebyeHuckel")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Debye-Huckel.");
        using real_type     = typename traitsT::real_type;
        using potential_t   = DebyeHuckelPotential<real_type>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        auto pot_para = read_debye_huckel_potential<traitsT>(global);

        return make_unique<interaction_t>(
                std::move(pot_para.first), std::move(pot_para.second),
                read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "LennardJones")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Lennard-Jones.");
        using real_type     = typename traitsT::real_type;
        using potential_t   = LennardJonesPotential<real_type>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        auto pot_para = read_lennard_jones_potential<traitsT>(global);

        return make_unique<interaction_t>(
                std::move(pot_para.first), std::move(pot_para.second),
                read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "UniformLennardJones")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is uniform Lennard-Jones.");
        using real_type     = typename traitsT::real_type;
        using potential_t   = UniformLennardJonesPotential<real_type>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        auto pot_para = read_uniform_lennard_jones_potential<traitsT>(global);

        return make_unique<interaction_t>(
                std::move(pot_para.first), std::move(pot_para.second),
                read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "LennardJonesAttractive")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is the attractive part of Lennard-Jones.");
        using real_type     = typename traitsT::real_type;
        using potential_t   = LennardJonesAttractivePotential<real_type>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        auto pot_para = read_lennard_jones_attractive_potential<traitsT>(global);

        return make_unique<interaction_t>(
                std::move(pot_para.first), std::move(pot_para.second),
                read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "WCA")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is WCA.");
        using real_type     = typename traitsT::real_type;
        using potential_t   = WCAPotential<real_type>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        auto pot_para = read_wca_potential<traitsT>(global);

        return make_unique<interaction_t>(
                std::move(pot_para.first), std::move(pot_para.second),
                read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "3SPN2ExcludedVolume")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is 3SPN2ExcludedVolume.");
        using real_type     = typename traitsT::real_type;
        using potential_t   = ThreeSPN2ExcludedVolumePotential<real_type>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        auto pot_para = read_3spn2_excluded_volume_potential<traitsT>(global);

        return make_unique<interaction_t>(
                std::move(pot_para.first), std::move(pot_para.second),
                read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "iSoLFAttractive")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is iSoLFAttractive.");
        using real_type     = typename traitsT::real_type;
        using potential_t   = iSoLFAttractivePotential<real_type>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        auto pot_para = read_isolf_potential<traitsT>(global);

        return make_unique<interaction_t>(
                std::move(pot_para.first), std::move(pot_para.second),
                read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "UniformCubicPan")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is UniformCubicPan");
        using real_type     = typename traitsT::real_type;
        using potential_t   = UniformCubicPanPotential<real_type>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        auto pot_para = read_uniform_cubic_pan_potential<traitsT>(global);

        return make_unique<interaction_t>(
            std::move(pot_para.first), std::move(pot_para.second),
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
            "- \"3SPN2ExcludedVolume\"  : excluded volume for 3SPN2 DNA model",
            "- \"iSoLFAttractive\"      : attractive potential for iSoLF lipid model",
            "- \"UniformCubicPan\" : cubic function potential with uniform paramters"
            }));
    }
}

// ----------------------------------------------------------------------------
// 3SPN2 Base-Base Interaction

template<typename traitsT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_3spn2_base_pair_interaction(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type      = typename traitsT::real_type;
    using base_kind      = parameter_3SPN2::base_kind;
    using parameter_list = ThreeSPN2BasePairParameterList<traitsT>;
    using potential_type = ThreeSPN2BasePairPotential<real_type>;
    using parameter_type = typename parameter_list::parameter_type;

    // [[forcefields.global]]
    // interaction = "3SPN2BasePair"
    // potential   = "3SPN2"
    // spatial_partition = {type = "CellList", margin = 1.0}
    // parameters = [
    // {strand = 0,          S =   0, B =   1, Base = "A"},
    // {strand = 0, P =   2, S =   3, B =   4, Base = "T"},
    // # ...
    // {strand = 0, P =  92, S =  93, B =  94, Base = "C"},
    // {strand = 1,          S =  95, B =  96, Base = "G"},
    // # ...
    // {strand = 1, P = 184, S = 185, B = 186, Base = "A"},
    // {strand = 1, P = 187, S = 188, B = 189, Base = "T"},
    // ]

    // ------------------------------------------------------------------------
    // read parameters

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    using nucleotide_info_type = parameter_3SPN2::NucleotideInfo;
    std::vector<nucleotide_info_type> nuc_idxs;

    nuc_idxs.reserve(ps.size());
    for(const auto& item : ps)
    {
        nucleotide_info_type nuc_idx;
        const auto ofs = find_parameter_or<std::int64_t>(item, env, "offset", 0);

        // at the edge of the DNA, Phosphate may not exist.
        if(item.as_table().count("P") != 0)
        {
            nuc_idx.P = find_parameter<std::size_t>(item, env, "P") + ofs;
        }
        nuc_idx.S          = find_parameter<std::size_t>(item, env, "S") + ofs;
        nuc_idx.B          = find_parameter<std::size_t>(item, env, "B") + ofs;
        nuc_idx.strand     = find_parameter<std::size_t>(item, env, "strand");

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

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(nuc_idxs.size());
    for(std::size_t i=0; i<nuc_idxs.size(); ++i)
    {
        const auto& base = nuc_idxs.at(i);
        const auto  B = base.B;

        parameter_type p;
        p.base   = base.base;
        p.S_idx  = base.S;

        MJOLNIR_LOG_INFO("Base idx = ", B, ", base = ", p.base, ", Sugar idx = ", p.S_idx);

        params.emplace_back(B, p);
    }

    const auto pot = toml::find<std::string>(global, "potential");
    if(pot == "3SPN2")
    {
        return make_unique<ThreeSPN2BasePairInteraction<traitsT>>(
            potential_type{},
            parameter_list(ThreeSPN2BasePairGlobalPotentialParameter<real_type>{},
                std::move(params), read_ignore_particles_within(global),
                read_ignored_molecule(global), read_ignored_group(global)),
            read_spatial_partition<traitsT, potential_type>(global));
    }
    else if(pot == "3SPN2C")
    {
        return make_unique<ThreeSPN2BasePairInteraction<traitsT>>(
            potential_type{},
            parameter_list(ThreeSPN2CBasePairGlobalPotentialParameter<real_type>{},
                std::move(params), read_ignore_particles_within(global),
                read_ignored_molecule(global), read_ignored_group(global)),
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

template<typename traitsT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_3spn2_cross_stacking_interaction(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type      = typename traitsT::real_type;
    using base_kind      = parameter_3SPN2::base_kind;
    using parameter_list = ThreeSPN2CrossStackingParameterList<traitsT>;
    using potential_type = ThreeSPN2CrossStackingPotential<real_type>;
    using parameter_type = typename parameter_list::parameter_type;

    // [[forcefields.global]]
    // interaction = "3SPN2CrossStacking"
    // potential   = "3SPN2"
    // spatial_partition = {type = "CellList", margin = 1.0}
    // parameters = [
    // {strand = 0, nucleotide =  0,          S =   0, B =   1, Base = "A"},
    // {strand = 0, nucleotide =  1, P =   2, S =   3, B =   4, Base = "T"},
    // # ...
    // {strand = 0, nucleotide = 31, P =  92, S =  93, B =  94, Base = "C"},
    // {strand = 1, nucleotide = 32,          S =  95, B =  96, Base = "G"},
    // # ...
    // {strand = 1, nucleotide = 62, P = 184, S = 185, B = 186, Base = "A"},
    // {strand = 1, nucleotide = 63, P = 187, S = 188, B = 189, Base = "T"},
    // ]

    // ------------------------------------------------------------------------
    // read parameters

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    using nucleotide_info_type = parameter_3SPN2::NucleotideInfo;
    std::vector<nucleotide_info_type> nuc_idxs;

    nuc_idxs.reserve(ps.size());
    for(const auto& item : ps)
    {
        nucleotide_info_type nuc_idx;
        const auto ofs = find_parameter_or<std::int64_t>(item, env, "offset", 0);

        // at the edge of the DNA, Phosphate may not exist.
        if(item.as_table().count("P") != 0)
        {
            nuc_idx.P = find_parameter<std::size_t>(item, env, "P") + ofs;
        }
        nuc_idx.S          = find_parameter<std::size_t>(item, env, "S") + ofs;
        nuc_idx.B          = find_parameter<std::size_t>(item, env, "B") + ofs;
        nuc_idx.strand     = find_parameter<std::size_t>(item, env, "strand");

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

        MJOLNIR_LOG_INFO("ThreeSPN2CrossStacking: nucleotide = ", nuc_idx);
        nuc_idxs.push_back(nuc_idx);
    }

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(nuc_idxs.size());
    for(std::size_t i=0; i<nuc_idxs.size(); ++i)
    {
        const auto& base = nuc_idxs.at(i);
        const auto  B = base.B;

        parameter_type p;
        p.base = base.base;
        p.S    = base.S;
        p.B3   = potential_type::invalid();
        p.B5   = potential_type::invalid();

        if (i+1 < nuc_idxs.size() &&
            nuc_idxs.at(i+1).strand == base.strand)
        {
            p.B3 = nuc_idxs.at(i+1).B;
        }
        if (0 < i && nuc_idxs.at(i-1).strand == base.strand)
        {
            p.B5 = nuc_idxs.at(i-1).B;
        }
        MJOLNIR_LOG_INFO("Base idx = ", B, ", base = ", p.base,
                ", Sugar idx = ", p.S, ", B5_idx = ", p.B5, ", B3_idx = ", p.B3);

        params.emplace_back(B, p);
    }

    const auto pot = toml::find<std::string>(global, "potential");
    if(pot == "3SPN2")
    {
        return make_unique<ThreeSPN2CrossStackingInteraction<traitsT>>(
            potential_type{},
            parameter_list(ThreeSPN2CrossStackingGlobalPotentialParameter<real_type>{},
                std::move(params), read_ignore_particles_within(global),
                read_ignored_molecule(global), read_ignored_group(global)),
            read_spatial_partition<traitsT, potential_type>(global));
    }
    else if(pot == "3SPN2C")
    {
        return make_unique<ThreeSPN2CrossStackingInteraction<traitsT>>(
            potential_type{},
            parameter_list(ThreeSPN2CCrossStackingGlobalPotentialParameter<real_type>{},
                std::move(params), read_ignore_particles_within(global),
                read_ignored_molecule(global), read_ignored_group(global)),
            read_spatial_partition<traitsT, potential_type>(global));
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_local_3spn2_cross_stacking_interaction: "
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
    using real_type           = typename traitsT::real_type;
    using parameter_list_type = ProteinDNANonSpecificParameterList<traitsT>;
    using contact_parameter_type = typename parameter_list_type::contact_parameter_type;
    using dna_index_type         = typename parameter_list_type::dna_index_type;
    using potential_type = ProteinDNANonSpecificPotential<real_type>;

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

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

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
    return make_unique<ProteinDNANonSpecificInteraction<traitsT>>(potential_type{},
        parameter_list_type(sgm, dlt, cutoff, std::move(contacts), std::move(dnas),
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
    using real_type              = typename traitsT::real_type;
    using parameter_list_type    = PWMcosParameterList<traitsT>;
    using contact_parameter_type = typename parameter_list_type::contact_parameter_type;
    using dna_parameter_type     = typename parameter_list_type::dna_parameter_type;
    using base_kind              = typename parameter_list_type::base_kind;
    using potential_type         = PWMcosPotential<real_type>;

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

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

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
    return make_unique<PWMcosInteraction<traitsT>>(potential_type{},
        parameter_list_type(sgm, phi, Eunit, Eshift, cutoff,
            std::move(contacts), std::move(dnas),
            read_ignore_particles_within(global), read_ignored_molecule(global),
            read_ignored_group(global)),
        read_spatial_partition<traitsT, potential_type>(global));
}

// ----------------------------------------------------------------------------
// Stoichiometric Interaction

template<typename traitsT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_stoichiometric_interaction(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    using real_type           = typename traitsT::real_type;
    using potential_type      = GlobalStoichiometricInteractionPotential<real_type>;
    using parameter_list_type = StoichiometricInteractionRule<traitsT, potential_type>;

    // ```toml
    // [[forcefields.global]]
    // interaction = "Stoichiometric"
    // potential   = "Stoichiometric"
    // epsilon     = 1.0 # the depth of potential valley
    // coefA = 3
    // coefB = 1
    // v0 = 4.0
    // parameters = [
    // {index = 0, kind = "A"},
    // {index = 1, kind = "B"},
    // # ...
    // ]
    // ```

    const auto        pot     = toml::find<std::string>(global, "potential");
    const real_type   epsilon = toml::find<real_type>  (global, "epsilon");
    const std::size_t coefa   = toml::find<std::size_t>(global, "coefA");
    const std::size_t coefb   = toml::find<std::size_t>(global, "coefB");
    const real_type   v0      = toml::find<real_type>  (global, "v0");
    const real_type   range   = toml::find<real_type>  (global, "range");
    const auto&       ps      = toml::find<toml::array>(global, "parameters");

    std::vector<std::size_t> a_indices;
    std::vector<std::size_t> b_indices;

    for(const auto& item : ps)
    {
        const std::size_t idx  = toml::find<std::size_t>(item, "index");
        const std::string kind = toml::find<std::string>(item, "kind");

        if(kind == "A")
        {
            a_indices.push_back(idx);
        }
        else if(kind == "B")
        {
            b_indices.push_back(idx);
        }
        else
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_stoichiometric_interaction: unknown particle kind ",
                toml::find<toml::value>(global, "kinds"), "here", {
                "expected value is \"A\" or \"B\""
                }));
        }
    }

    parameter_list_type parameter_list(
            std::move(a_indices), std::move(b_indices),
            read_ignore_particles_within(global),
            read_ignored_molecule(global), read_ignored_group(global));

    return make_unique<GlobalStoichiometricInteraction<traitsT, potential_type>>(
            potential_type{v0, range}, std::move(parameter_list),
            read_spatial_partition<traitsT, potential_type>(global),
            epsilon, coefa, coefb);
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
    else if(interaction == "3SPN2BasePair")
    {
        MJOLNIR_LOG_NOTICE("3SPN2BasePairInteraction found.");
        return read_global_3spn2_base_pair_interaction<traitsT>(global);
    }
    else if(interaction == "3SPN2CrossStacking")
    {
        MJOLNIR_LOG_NOTICE("3SPN2CrossStackingBasePairInteraction found.");
        return read_global_3spn2_cross_stacking_interaction<traitsT>(global);
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
    else if(interaction == "Stoichiometric")
    {
        MJOLNIR_LOG_NOTICE("Stoichiometric Interaction found.");
        return read_global_stoichiometric_interaction<traitsT>(global);
    }
    else
    {
        throw std::runtime_error(toml::format_error("[error] "
            "mjolnir::read_global_interaction: invalid interaction",
            toml::find<toml::value>(global, "interaction"), "here", {
            "expected value is one of the following.",
            "- \"Pair\"         : well-known pair interaction depends only on the distance",
            "- \"3SPN2BaseBase\": Base pair and cross stacking interaction for 3SPN2 DNA model",
            "- \"PDNS\"         : direction-dependent contact to represent hydrogen bond",
            "- \"PWMcos\"       : direction-dependent contact to reproduce DNA sequeucne-specific interaction"
            }));
    }
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_global_interaction(const toml::value&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_global_interaction(const toml::value&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_global_interaction(const toml::value&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_global_interaction(const toml::value&);

extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_global_3spn2_base_pair_interaction(const toml::value&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_global_3spn2_base_pair_interaction(const toml::value&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_global_3spn2_base_pair_interaction(const toml::value&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_global_3spn2_base_pair_interaction(const toml::value&);

extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_global_3spn2_cross_stacking_interaction(const toml::value&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_global_3spn2_cross_stacking_interaction(const toml::value&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_global_3spn2_cross_stacking_interaction(const toml::value&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_global_3spn2_cross_stacking_interaction(const toml::value&);

extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_global_pair_interaction(const toml::value&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_global_pair_interaction(const toml::value&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_global_pair_interaction(const toml::value&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_global_pair_interaction(const toml::value&);

extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_global_stoichiometric_interaction(const toml::value&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float , UnlimitedBoundary>       >> read_global_stoichiometric_interaction(const toml::value&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_global_stoichiometric_interaction(const toml::value&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float , CuboidalPeriodicBoundary>>> read_global_stoichiometric_interaction(const toml::value&);
#endif

} // mjolnir
#endif// MJOLNIR_READ_GLOBAL_INTERACTION_HPP
