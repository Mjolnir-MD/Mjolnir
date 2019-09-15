#ifndef MJOLNIR_INPUT_READ_GLOBAL_INTERACTION_HPP
#define MJOLNIR_INPUT_READ_GLOBAL_INTERACTION_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/interaction/global/GlobalPairInteraction.hpp>
#include <mjolnir/interaction/global/GlobalPairLennardJonesInteraction.hpp>
#include <mjolnir/interaction/global/GlobalPairUniformLennardJonesInteraction.hpp>
#include <mjolnir/interaction/global/GlobalPairExcludedVolumeInteraction.hpp>
#include <mjolnir/forcefield/3SPN2/ThreeSPN2BaseBaseInteraction.hpp>
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
    using real_type = typename traitsT::real_type;

    const auto potential = toml::find<std::string>(global, "potential");

    if(potential == "ExcludedVolume")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Excluded Volume.");
        using potential_t   = ExcludedVolumePotential<real_type>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        return make_unique<interaction_t>(
            read_excluded_volume_potential<real_type>(global),
            read_spatial_partition<traitsT, potential_t>(global));
    }
    if(potential == "HardCoreExcludedVolume")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Hard Core Excluded Volume.");
        using potential_t   = HardCoreExcludedVolumePotential<real_type>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        return make_unique<interaction_t>(
            read_hard_core_excluded_volume_potential<real_type>(global),
            read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "DebyeHuckel")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Debye-Huckel.");
        using potential_t   = DebyeHuckelPotential<real_type>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        return make_unique<interaction_t>(
            read_debye_huckel_potential<real_type>(global),
            read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "LennardJones")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Lennard-Jones.");
        using potential_t   = LennardJonesPotential<real_type>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        return make_unique<interaction_t>(
            read_lennard_jones_potential<real_type>(global),
            read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "UniformLennardJones")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Uniform Lennard-Jones.");
        using potential_t   = UniformLennardJonesPotential<real_type>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        return make_unique<interaction_t>(
            read_uniform_lennard_jones_potential<real_type>(global),
            read_spatial_partition<traitsT, potential_t>(global));
    }
    else if(potential == "3SPN2ExcludedVolume")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is 3SPN2ExcludedVolume.");
        using potential_t   = ThreeSPN2ExcludedVolumePotential<real_type>;
        using interaction_t = GlobalPairInteraction<traitsT, potential_t>;

        return make_unique<interaction_t>(
            read_3spn2_excluded_volume_potential<real_type>(global),
            read_spatial_partition<traitsT, potential_t>(global));
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_global_pair_interaction: invalid potential",
            toml::find<toml::value>(global, "potential"), "here", {
            "expected value is one of the following.",
            "- \"ExcludedVolume\"       : repulsive r^12 potential",
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
    using potential_type      = ThreeSPN2BaseBaseInteractionPotential<real_type>;
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

        // at the edge of the DNA, Phosphate may not exist.
        if(item.as_table().count("P") != 0)
        {
            nuc_idx.P = find_parameter<std::size_t>(item, env, "P");
        }
        nuc_idx.S          = find_parameter<std::size_t>(item, env, "S");
        nuc_idx.B          = find_parameter<std::size_t>(item, env, "B");
        nuc_idx.strand     = find_parameter<std::size_t>(item, env, "strand");
        nuc_idx.nucleotide = find_parameter<std::size_t>(item, env, "nucleotide");

        const auto bk      = find_parameter<std::string>(item, env, "Base");
        if     (bk == "A") {nuc_idx.base = base_kind::A;}
        else if(bk == "T") {nuc_idx.base = base_kind::T;}
        else if(bk == "G") {nuc_idx.base = base_kind::G;}
        else if(bk == "C") {nuc_idx.base = base_kind::C;}
        else
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_local_3spn2_base_stacking_interaction: "
                "invalid Base", item, "here", {
                "expected value is one of the \"A\", \"T\", \"C\", \"G\"."
                }));
        }

        MJOLNIR_LOG_INFO("ThreeSPN2BaseStackingPotential: nucleotide = ", nuc_idx);
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
    if(interaction == "3SPN2BaseBase")
    {
        MJOLNIR_LOG_NOTICE("3SPN2BaseBaseInteraction found.");
        return read_global_3spn2_base_base_interaction<traitsT>(global);
    }
    else
    {
        throw std::runtime_error(toml::format_error("[error] "
            "mjolnir::read_global_interaction: invalid interaction",
            toml::find<toml::value>(global, "interaction"), "here", {
            "expected value is one of the following.",
            "- \"Pair\": well-known pair interaction depends only on the distance",
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
