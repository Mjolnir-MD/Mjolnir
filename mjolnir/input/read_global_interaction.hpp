#ifndef MJOLNIR_READ_GLOBAL_INTERACTION_HPP
#define MJOLNIR_READ_GLOBAL_INTERACTION_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/interaction/GlobalPairInteraction.hpp>
#include <mjolnir/interaction/specialization/GlobalPairLennardJonesInteraction.hpp>
#include <mjolnir/interaction/specialization/GlobalPairUniformLennardJonesInteraction.hpp>
#include <mjolnir/interaction/specialization/GlobalPairExcludedVolumeInteraction.hpp>
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
    MJOLNIR_SCOPE(read_global_pair_interaction(), 0);
    using real_type = typename traitsT::real_type;

    const auto potential = toml::find<std::string>(global, "potential");

    if(potential == "ExcludedVolume")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Excluded Volume.");
        using potential_t = ExcludedVolumePotential<real_type>;

        return read_spatial_partition<traitsT, potential_t>(
            global, read_excluded_volume_potential<real_type>(global));
    }
    else if(potential == "DebyeHuckel")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Debye-Huckel.");
        using potential_t = DebyeHuckelPotential<real_type>;

        return read_spatial_partition<traitsT, potential_t>(
            global, read_debye_huckel_potential<real_type>(global));
    }
    else if(potential == "LennardJones")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Lennard-Jones.");
        using potential_t = LennardJonesPotential<real_type>;

        return read_spatial_partition<traitsT, potential_t>(
            global, read_lennard_jones_potential<real_type>(global));
    }
    else if(potential == "UniformLennardJones")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Uniform Lennard-Jones.");
        using potential_t = UniformLennardJonesPotential<real_type>;

        return read_spatial_partition<traitsT, potential_t>(
            global, read_uniform_lennard_jones_potential<real_type>(global));
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
            "- \"UniformLennardJones\"  : famous r^12 - r^6 potential with uniform parameters"
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
    MJOLNIR_SCOPE(read_global_interaction(), 0);
    const auto interaction = toml::find<std::string>(global, "interaction");

    if(interaction == "Pair")
    {
        MJOLNIR_LOG_NOTICE("Pair interaction found.");
        return read_global_pair_interaction<traitsT>(global);
    }
    else
    {
        throw std::runtime_error(toml::format_error("[error] "
            "mjolnir::read_global_interaction: invalid interaction",
            toml::find<toml::value>(global, "interaction"), "here", {
            "expected value is one of the following.",
            "- \"Pair\": well-known pair interaction depends only on the distance"
            }));
    }
}

} // mjolnir
#endif// MJOLNIR_READ_GLOBAL_INTERACTION_HPP
