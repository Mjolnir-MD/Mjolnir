#ifndef MJOLNIR_READ_EXTERNAL_INTERACTION_HPP
#define MJOLNIR_READ_EXTERNAL_INTERACTION_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/interaction/ExternalDistanceInteraction.hpp>
#include <mjolnir/core/AxisAlignedPlane.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/get_toml_value.hpp>
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
read_external_distance_interaction(const toml::Table& external, shapeT&& shape)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_global_distance_interaction(), 0);
    using real_type = typename traitsT::real_type;

    const auto potential = get_toml_value<std::string>(
            external, "potential", "[forcefield.external]");

    if(potential == "ImplicitMembrane")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Implicit Membrane.");
        using potential_t   = ImplicitMembranePotential<real_type>;
        using interaction_t = ExternalDistanceInteraction<
                                    traitsT, potential_t, shapeT>;

        return make_unique<interaction_t>(std::move(shape),
             read_implicit_membrane_potential<real_type>(external));
    }
    else if(potential == "LennardJonesWall")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Lennard-Jones.");
        using potential_t   = LennardJonesWallPotential<real_type>;
        using interaction_t = ExternalDistanceInteraction<
                                    traitsT, potential_t, shapeT>;
        return make_unique<interaction_t>(std::move(shape),
             read_lennard_jones_wall_potential<real_type>(external));
    }
    else if(potential == "ExcludedVolumeWall")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Excluded-Volume.");
        using potential_t   = ExcludedVolumeWallPotential<real_type>;
        using interaction_t = ExternalDistanceInteraction<
                                    traitsT, potential_t, shapeT>;
        return make_unique<interaction_t>(std::move(shape),
             read_excluded_volume_wall_potential<real_type>(external));
    }
    else
    {
        throw_exception<std::runtime_error>(
            "invalid potential as ExternalDistanceInteraction: ", potential);
    }
}

template<typename traitsT>
std::unique_ptr<ExternalForceInteractionBase<traitsT>>
read_external_distance_interaction_shape(const toml::Table& external)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_global_distance_interaction(), 0);
    using real_type = typename traitsT::real_type;

    const auto shape = get_toml_value<toml::Table>(external, "shape",
            "[forcefield.external] for ExternalDistance");
    const auto name  = get_toml_value<std::string>(shape, "name",
            "[forcefield.external.shape] for ExternalDistance");

    if(name == "AxisAlignedPlane")
    {
        const auto pos = get_toml_value<real_type>(shape, "position",
            "[forcefield.external.shape] for ExternalDistance");
        const auto margin = get_toml_value<real_type>(shape, "margin",
            "[forcefield.external.shape] for ExternalDistance");
        MJOLNIR_LOG_INFO("pos    = ", pos);
        MJOLNIR_LOG_INFO("margin = ", margin);

        const auto axis = get_toml_value<std::string>(shape, "axis",
            "[forcefield.external.shape] for ExternalDistance");
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
            throw std::runtime_error("invalid axis name: " + axis);
        }
    }
    else
    {
        throw std::runtime_error("invalid external forcefield shape: " + name);
    }
}

// ----------------------------------------------------------------------------
// general read_(local|global|external)_interaction function
// ----------------------------------------------------------------------------

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_local_interaction(const toml::Table& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_local_interaction(), 0);
    const auto interaction = get_toml_value<std::string>(
            local, "interaction", "[forcefields.local]");

    const auto kind = get_toml_value<std::string>(
            local, "topology", "[forcefield.local]");
    MJOLNIR_LOG_INFO("connection kind = ", kind);

    if(interaction == "BondLength")
    {
        MJOLNIR_LOG_NOTICE("Bond Length interaction found.");
        return read_bond_length_interaction<traitsT>(kind, local);
    }
    else if(interaction == "BondAngle")
    {
        MJOLNIR_LOG_NOTICE("Bond Angle interaction found.");
        return read_bond_angle_interaction<traitsT>(kind, local);
    }
    else if(interaction == "DihedralAngle")
    {
        MJOLNIR_LOG_NOTICE("Dihedral Angle interaction found.");
        return read_dihedral_angle_interaction<traitsT>(kind, local);
    }
    else
    {
        throw std::runtime_error(
                "invalid local interaction type: " + interaction);
    }
}

template<typename traitsT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_interaction(const toml::Table& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_global_interaction(), 0);
    const auto interaction = get_toml_value<std::string>(
            global, "interaction", "[forcefields.global]");

    if(interaction == "Pair")
    {
        MJOLNIR_LOG_NOTICE("Pair interaction found.");
        return read_global_pair_interaction<traitsT>(global);
    }
    else
    {
        throw std::runtime_error(
                "invalid global interaction type: " + interaction);
    }
}

template<typename traitsT>
std::unique_ptr<ExternalForceInteractionBase<traitsT>>
read_external_interaction(const toml::Table& external)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_external_interaction(), 0);
    const auto interaction = get_toml_value<std::string>(
            external, "interaction", "[forcefields.external]");

    if(interaction == "Distance")
    {
        MJOLNIR_LOG_NOTICE("Distance interaction found.");
        return read_external_distance_interaction_shape<traitsT>(external);
    }
    else
    {
        throw std::runtime_error(
                "invalid external interaction type: " + interaction);
    }
}

} // mjolnir
#endif// MJOLNIR_READ_EXTERNAL_INTERACTION_HPP
