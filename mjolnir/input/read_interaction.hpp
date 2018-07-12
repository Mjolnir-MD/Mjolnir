#ifndef MJOLNIR_READ_INTERACTION
#define MJOLNIR_READ_INTERACTION
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/BondLengthInteraction.hpp>
#include <mjolnir/core/BondAngleInteraction.hpp>
#include <mjolnir/core/DihedralAngleInteraction.hpp>
#include <mjolnir/core/GlobalDistanceInteraction.hpp>
#include <mjolnir/core/AxisAlignedPlane.hpp>
#include <mjolnir/core/ExternalDistanceInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_potential.hpp>
#include <mjolnir/input/read_spatial_partition.hpp>
#include <memory>

namespace mjolnir
{

// ----------------------------------------------------------------------------
// local interaction
// ----------------------------------------------------------------------------

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_bond_length_interaction(const std::string& kind, const toml::Table& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_bond_length_interaction(), 0);

    const auto potential = get_toml_value<std::string>(
            local, "potential", "[forcefield.local]");

    if(potential == "Harmonic")
    {
        MJOLNIR_SCOPE(potential == "Harmonic", 1);
        using potentialT = HarmonicPotential<traitsT>;

        return make_unique<BondLengthInteraction<traitsT, potentialT>>(
                kind, read_local_potential<potentialT, 2>(local));
    }
    else if(potential == "Go1012Contact")
    {
        MJOLNIR_SCOPE(potential == "Go1012Contact", 1);
        using potentialT = Go1012ContactPotential<traitsT>;

        return make_unique<BondLengthInteraction<traitsT, potentialT>>(
                kind, read_local_potential<potentialT, 2>(local));
    }
    else if(potential == "Gaussian")
    {
        MJOLNIR_SCOPE(potential == "Gaussian", 1);
        using potentialT = GaussianPotential<traitsT>;

        return make_unique<BondLengthInteraction<traitsT, potentialT>>(
                kind, read_local_potential<potentialT, 2>(local));
    }
    else
    {
        throw_exception<std::runtime_error>(
                "invalid potential as BondLengthInteraction: ", potential);
    }
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_bond_angle_interaction(
    const typename LocalInteractionBase<traitsT>::connection_kind_type kind,
    const toml::Table& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_bond_angle_interaction(), 0);

    const auto potential = get_toml_value<std::string>(
            local, "potential", "[[forcefields.local]]");

    if(potential == "Harmonic")
    {
        MJOLNIR_SCOPE(potential == "Harmonic", 1);
        using potentialT = HarmonicPotential<traitsT>;

        return make_unique<BondAngleInteraction<traitsT, potentialT>>(
                kind, read_local_potential<potentialT, 3>(local));
    }
    else if(potential == "FlexibleLocalAngle")
    {
        MJOLNIR_SCOPE(potential == "FlexibleLocalAngle", 1);
        using potentialT = FlexibleLocalAnglePotential<traitsT>;

        return make_unique<BondAngleInteraction<traitsT, potentialT>>(
                kind, read_local_potential<potentialT, 3>(local));
    }
    else if(potential == "Gaussian")
    {
        MJOLNIR_SCOPE(potential == "Gaussian", 1);
        using potentialT = GaussianPotential<traitsT>;

        return make_unique<BondAngleInteraction<traitsT, potentialT>>(
                kind, read_local_potential<potentialT, 3>(local));
    }
    else
    {
        throw_exception<std::runtime_error>(
                "invalid potential as BondAngleInteraction: " + potential);
    }
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_dihedral_angle_interaction(
    const typename LocalInteractionBase<traitsT>::connection_kind_type kind,
    const toml::Table& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_dihedral_angle_interaction(), 0);

    const auto potential = get_toml_value<std::string>(
            local, "potential", "[[forcefields.local]]");

    if(potential == "Harmonic")
    {
        MJOLNIR_SCOPE(potential == "Harmonic", 1);
        using potentialT = HarmonicPotential<traitsT>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<potentialT, 4>(local));
    }
    else if(potential == "ClementiDihedral")
    {
        MJOLNIR_SCOPE(potential == "ClementiDihedral", 1);
        using potentialT = ClementiDihedralPotential<traitsT>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<potentialT, 4>(local));
    }
    else if(potential == "Gaussian")
    {
        MJOLNIR_SCOPE(potential == "Gaussian", 1);
        using potentialT = GaussianPotential<traitsT>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<potentialT, 4>(local));
    }
    else if(potential == "FlexibleLocalDihedral")
    {
        MJOLNIR_SCOPE(potential == "FlexibleLocalDihedral", 1);
        using potentialT = FlexibleLocalDihedralPotential<traitsT>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<potentialT, 4>(local));
    }
    else
    {
        throw_exception<std::runtime_error>(
                "invalid potential as DihedralAngleInteraction: " + potential);
    }
}

// ----------------------------------------------------------------------------
// global interaction
// ----------------------------------------------------------------------------

template<typename traitsT, typename ignoreT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_distance_interaction(const toml::Table& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_global_distance_interaction(), 0);

    const auto potential = get_toml_value<std::string>(
            global, "potential", "[[forcefield.global]]");

    if(potential == "ExcludedVolume")
    {
        MJOLNIR_SCOPE(potential == "ExcludedVolume", 1);
        using potential_t = ExcludedVolumePotential<traitsT, ignoreT>;

        return read_spatial_partition_for_distance<traitsT, potential_t>(
            global, read_excluded_volume_potential<traitsT, ignoreT>(global));
    }
    else if(potential == "DebyeHuckel")
    {
        MJOLNIR_SCOPE(potential == "DebyeHuckel", 1);
        using potential_t = DebyeHuckelPotential<traitsT, ignoreT>;

        return read_spatial_partition_for_distance<traitsT, potential_t>(
            global, read_debye_huckel_potential<traitsT, ignoreT>(global));
    }
    else if(potential == "LennardJones")
    {
        MJOLNIR_SCOPE(potential == "LennardJones", 1);
        using potential_t = LennardJonesPotential<traitsT, ignoreT>;

        return read_spatial_partition_for_distance<traitsT, potential_t>(
            global, read_lennard_jones_potential<traitsT, ignoreT>(global));
    }
    else
    {
        throw_exception<std::runtime_error>(
                "invalid potential as GlobalDistanceInteraction: ", potential);
    }
}

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
        MJOLNIR_SCOPE(potential == "ImplicitMembrane", 1);
        using potential_t   = ImplicitMembranePotential<traitsT>;
        using interaction_t = ExternalDistanceInteraction<
                                    traitsT, potential_t, shapeT>;

        return make_unique<interaction_t>(std::move(shape),
             read_implicit_membrane_potential<traitsT>(external));
    }
    else if(potential == "LennardJonesWall")
    {
        MJOLNIR_SCOPE(potential == "LennardJonesWall", 1);
        using potential_t   = LennardJonesWallPotential<traitsT>;
        using interaction_t = ExternalDistanceInteraction<
                                    traitsT, potential_t, shapeT>;
        return make_unique<interaction_t>(std::move(shape),
             read_lennard_jones_wall_potential<traitsT>(external));
    }
    else if(potential == "ExcludedVolumeWall")
    {
        MJOLNIR_SCOPE(potential == "ExcludedVolumeWall", 1);
        using potential_t   = ExcludedVolumeWallPotential<traitsT>;
        using interaction_t = ExternalDistanceInteraction<
                                    traitsT, potential_t, shapeT>;
        return make_unique<interaction_t>(std::move(shape),
             read_excluded_volume_wall_potential<traitsT>(external));
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
        MJOLNIR_SCOPE(potential == "AxisAlignedPlane", 1);
        const auto pos = get_toml_value<real_type>(shape, "position",
            "[forcefield.external.shape] for ExternalDistance");
        const auto margin = get_toml_value<real_type>(shape, "margin",
            "[forcefield.external.shape] for ExternalDistance");
        MJOLNIR_LOG_INFO("pos    = ", pos);
        MJOLNIR_LOG_INFO("margin = ", margin);

        const auto axis = get_toml_value<std::string>(shape, "axis",
            "[forcefield.external.shape] for ExternalDistance");
        if(axis == "X" || axis == "+X")
        {
            MJOLNIR_LOG_INFO("axis   = +X");
            using shape_t = AxisAlignedPlane<traitsT, PositiveXDirection>;
            return read_external_distance_interaction<traitsT, shape_t>(
                    external, shape_t(pos, margin));
        }
        else if(axis == "-X")
        {
            MJOLNIR_LOG_INFO("axis   = -X");
            using shape_t = AxisAlignedPlane<traitsT, NegativeXDirection>;
            return read_external_distance_interaction<traitsT, shape_t>(
                    external, shape_t(pos, margin));
        }
        else if(axis == "Y" || axis == "+Y")
        {
            MJOLNIR_LOG_INFO("axis   = +Y");
            using shape_t = AxisAlignedPlane<traitsT, PositiveYDirection>;
            return read_external_distance_interaction<traitsT, shape_t>(
                    external, shape_t(pos, margin));
        }
        else if(axis == "-Y")
        {
            MJOLNIR_LOG_INFO("axis   = -Y");
            using shape_t = AxisAlignedPlane<traitsT, NegativeYDirection>;
            return read_external_distance_interaction<traitsT, shape_t>(
                    external, shape_t(pos, margin));
        }
        else if(axis == "Z" || axis == "+Z")
        {
            MJOLNIR_LOG_INFO("axis   = +Z");
            using shape_t = AxisAlignedPlane<traitsT, PositiveZDirection>;
            return read_external_distance_interaction<traitsT, shape_t>(
                    external, shape_t(pos, margin));
        }
        else if(axis == "-Z")
        {
            MJOLNIR_LOG_INFO("axis   = -Z");
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
        return read_bond_length_interaction<traitsT>(kind, local);
    }
    else if(interaction == "BondAngle")
    {
        return read_bond_angle_interaction<traitsT>(kind, local);
    }
    else if(interaction == "DihedralAngle")
    {
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
    const auto ignored_chain = get_toml_value<std::string>(
            global, "ignored_chain", "[forcefields.global]");

    if(interaction == "Distance")
    {
        MJOLNIR_LOG_INFO("Distance interaction found");
        if(ignored_chain == "Nothing")
        {
            MJOLNIR_LOG_INFO("all the interactions(both (inter|intra)-chain) are included");
            return read_global_distance_interaction<traitsT, IgnoreNothing>(global);
        }
        else if(ignored_chain == "Self")
        {
            MJOLNIR_LOG_INFO("intra-chain interaction is ignored");
            return read_global_distance_interaction<traitsT, IgnoreSelf>(global);
        }
        else if(ignored_chain == "Others")
        {
            MJOLNIR_LOG_INFO("inter-chain interaction is ignored");
            return read_global_distance_interaction<traitsT, IgnoreOthers>(global);
        }
        else
        {
            throw std::runtime_error("invalid `ignored_chain`: " + ignored_chain);
        }
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
        MJOLNIR_LOG_INFO("Distance interaction found");
        return read_external_distance_interaction_shape<traitsT>(external);
    }
    else
    {
        throw std::runtime_error(
                "invalid external interaction type: " + interaction);
    }
}

} // mjolnir
#endif// MJOLNIR_READ_INTERACTION
