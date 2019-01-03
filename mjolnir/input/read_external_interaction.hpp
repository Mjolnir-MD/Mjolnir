#ifndef MJOLNIR_READ_EXTERNAL_INTERACTION_HPP
#define MJOLNIR_READ_EXTERNAL_INTERACTION_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/interaction/ExternalDistanceInteraction.hpp>
#include <mjolnir/core/AxisAlignedPlane.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_external_potential.hpp>
#include <memory>

namespace mjolnir
{

// ----------------------------------------------------------------------------
// external interaction
// ----------------------------------------------------------------------------

template<typename traitsT, typename shapeT>
std::unique_ptr<ExternalForceInteractionBase<traitsT>>
read_external_distance_interaction(const toml::value& external, shapeT&& shape)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_global_distance_interaction(), 0);
    using real_type = typename traitsT::real_type;

    const auto potential = toml::find<std::string>(external, "potential");
    if(potential == "ImplicitMembrane")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Implicit Membrane.");
        using potential_t   = ImplicitMembranePotential<real_type>;
        using interaction_t = ExternalDistanceInteraction<
                                    traitsT, potential_t, shapeT>;

        return make_unique<interaction_t>(std::move(shape),
             read_implicit_membrane_potential<real_type>(external));
    }
    else if(potential == "LennardJonesWall")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Lennard-Jones.");
        using potential_t   = LennardJonesWallPotential<real_type>;
        using interaction_t = ExternalDistanceInteraction<
                                    traitsT, potential_t, shapeT>;
        return make_unique<interaction_t>(std::move(shape),
             read_lennard_jones_wall_potential<real_type>(external));
    }
    else if(potential == "ExcludedVolumeWall")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Excluded-Volume.");
        using potential_t   = ExcludedVolumeWallPotential<real_type>;
        using interaction_t = ExternalDistanceInteraction<
                                    traitsT, potential_t, shapeT>;
        return make_unique<interaction_t>(std::move(shape),
             read_excluded_volume_wall_potential<real_type>(external));
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_external_distance_interaction: invalid potential",
            toml::find<toml::value>(external, "potential"), "here", {
            "expected value is one of the following.",
            "- \"ExcludedVolumeWall\": repulsive r^12 potential",
            "- \"LennardJonesWall\"  : famous r^12 - r^6 potential",
            "- \"ImplicitMembrane\"  : single-well interaction that models a membrane"
            }));
    }
}

template<typename traitsT>
std::unique_ptr<ExternalForceInteractionBase<traitsT>>
read_external_distance_interaction_shape(const toml::value& external)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_global_distance_interaction(), 0);
    using real_type = typename traitsT::real_type;

    const auto shape = toml::find<toml::value>(external, "shape");
    const auto name  = toml::find<std::string>(shape, "name");

    if(name == "AxisAlignedPlane")
    {
        const auto pos    = toml::find<real_type>(shape, "position");
        const auto margin = toml::find<real_type>(shape, "margin");
        MJOLNIR_LOG_INFO("pos    = ", pos);
        MJOLNIR_LOG_INFO("margin = ", margin);

        const auto axis = toml::find<std::string>(shape, "axis");
        MJOLNIR_LOG_INFO("axis   = ", axis);

        if(axis == "X" || axis == "+X")
        {
            MJOLNIR_LOG_NOTICE("-- interaction shape is YZ-Plane (+).");
            using shape_t = AxisAlignedPlane<traitsT, PositiveXDirection>;
            return read_external_distance_interaction<traitsT, shape_t>(
                    external, shape_t(pos, margin));
        }
        else if(axis == "-X")
        {
            MJOLNIR_LOG_NOTICE("-- interaction shape is YZ-Plane (-).");
            using shape_t = AxisAlignedPlane<traitsT, NegativeXDirection>;
            return read_external_distance_interaction<traitsT, shape_t>(
                    external, shape_t(pos, margin));
        }
        else if(axis == "Y" || axis == "+Y")
        {
            MJOLNIR_LOG_NOTICE("-- interaction shape is ZX-Plane (+).");
            using shape_t = AxisAlignedPlane<traitsT, PositiveYDirection>;
            return read_external_distance_interaction<traitsT, shape_t>(
                    external, shape_t(pos, margin));
        }
        else if(axis == "-Y")
        {
            MJOLNIR_LOG_NOTICE("-- interaction shape is ZX-Plane (-).");
            using shape_t = AxisAlignedPlane<traitsT, NegativeYDirection>;
            return read_external_distance_interaction<traitsT, shape_t>(
                    external, shape_t(pos, margin));
        }
        else if(axis == "Z" || axis == "+Z")
        {
            MJOLNIR_LOG_NOTICE("-- interaction shape is XY-Plane (+).");
            using shape_t = AxisAlignedPlane<traitsT, PositiveZDirection>;
            return read_external_distance_interaction<traitsT, shape_t>(
                    external, shape_t(pos, margin));
        }
        else if(axis == "-Z")
        {
            MJOLNIR_LOG_NOTICE("-- interaction shape is XY-Plane (-).");
            using shape_t = AxisAlignedPlane<traitsT, NegativeZDirection>;
            return read_external_distance_interaction<traitsT, shape_t>(
                    external, shape_t(pos, margin));
        }
        else
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_external_distance_interaction_shape: invalid shape",
                toml::find<toml::value>(external, "shape"), "expected [+-]?[XYZ]"));
        }
    }
    else
    {
        throw std::runtime_error("invalid external forcefield shape: " + name);
    }
}

// ----------------------------------------------------------------------------
// general read_external_interaction function
// ----------------------------------------------------------------------------

template<typename traitsT>
std::unique_ptr<ExternalForceInteractionBase<traitsT>>
read_external_interaction(const toml::value& external)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_external_interaction(), 0);
    const auto interaction = toml::find<std::string>(external, "interaction");

    if(interaction == "Distance")
    {
        MJOLNIR_LOG_NOTICE("Distance interaction found.");
        return read_external_distance_interaction_shape<traitsT>(external);
    }
    else
    {
        throw std::runtime_error(toml::format_error("[error] "
            "mjolnir::read_external_interaction: invalid interaction",
            toml::find<toml::value>(external, "interaction"), "here", {
            "expected value is one of the following.",
            "- \"Distance\": interaction depending on the distance between particle and spatial structure"
            }));
    }
}

} // mjolnir
#endif// MJOLNIR_READ_EXTERNAL_INTERACTION_HPP
