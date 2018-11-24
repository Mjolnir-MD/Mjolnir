#ifndef MJOLNIR_READ_INTERACTION
#define MJOLNIR_READ_INTERACTION
#include <extlib/toml/toml.hpp>
#include <mjolnir/interaction/BondLengthInteraction.hpp>
#include <mjolnir/interaction/BondAngleInteraction.hpp>
#include <mjolnir/interaction/DihedralAngleInteraction.hpp>
#include <mjolnir/interaction/GlobalPairInteraction.hpp>
#include <mjolnir/interaction/ExternalDistanceInteraction.hpp>
#include <mjolnir/interaction/specialization/GlobalPairLennardJonesInteraction.hpp>
#include <mjolnir/interaction/specialization/GlobalPairUniformLennardJonesInteraction.hpp>
#include <mjolnir/interaction/specialization/GlobalPairExcludedVolumeInteraction.hpp>
#include <mjolnir/core/AxisAlignedPlane.hpp>
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
    using real_type = typename traitsT::real_type;

    const auto potential = get_toml_value<std::string>(
            local, "potential", "[forcefield.local]");

    if(potential == "Harmonic")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Harmonic.");
        using potentialT = HarmonicPotential<real_type>;

        return make_unique<BondLengthInteraction<traitsT, potentialT>>(
                kind, read_local_potential<2, potentialT>(local));
    }
    else if(potential == "Go1012Contact")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is 10-12 Go contact.");
        using potentialT = Go1012ContactPotential<real_type>;

        return make_unique<BondLengthInteraction<traitsT, potentialT>>(
                kind, read_local_potential<2, potentialT>(local));
    }
    else if(potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Gaussian.");
        using potentialT = GaussianPotential<real_type>;

        return make_unique<BondLengthInteraction<traitsT, potentialT>>(
                kind, read_local_potential<2, potentialT>(local));
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
    using real_type = typename traitsT::real_type;

    const auto potential = get_toml_value<std::string>(
            local, "potential", "[[forcefields.local]]");

    if(potential == "Harmonic")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Harmonic.");
        using potentialT = HarmonicPotential<real_type>;

        return make_unique<BondAngleInteraction<traitsT, potentialT>>(
                kind, read_local_potential<3, potentialT>(local));
    }
    else if(potential == "FlexibleLocalAngle")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Flexible Local Angle.");
        using potentialT = FlexibleLocalAnglePotential<real_type>;

        return make_unique<BondAngleInteraction<traitsT, potentialT>>(
                kind, read_local_potential<3, potentialT>(local));
    }
    else if(potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Gaussian.");
        using potentialT = GaussianPotential<real_type>;

        return make_unique<BondAngleInteraction<traitsT, potentialT>>(
                kind, read_local_potential<3, potentialT>(local));
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
    using real_type = typename traitsT::real_type;

    const auto potential = get_toml_value<std::string>(
            local, "potential", "[[forcefields.local]]");

    if(potential == "Harmonic")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Harmonic.");
        using potentialT = HarmonicPotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    else if(potential == "ClementiDihedral")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Clementi-Go's dihedral.");
        using potentialT = ClementiDihedralPotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    else if(potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Gaussian.");
        using potentialT = AngularGaussianPotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    else if(potential == "FlexibleLocalDihedral")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Flexible Local Dihedral.");
        using potentialT = FlexibleLocalDihedralPotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    // XXX generalization of this feature is too difficult (technically, it's
    // not so difficult, but practically, it makes the code messy...).
    else if(potential == "Gaussian+FlexibleLocalDihedral")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Gaussian + FlexibleLocalDihedral.");
        using potential_1_T = GaussianPotential<real_type>;
        using potential_2_T = FlexibleLocalDihedralPotential<real_type>;
        using potentialT    =
            SumLocalPotential<real_type, potential_1_T, potential_2_T>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(kind,
            read_local_potentials<4, real_type, potential_1_T, potential_2_T>(
                local, "Gaussian", "FlexibleLocalDihedral"));
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

template<typename traitsT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_pair_interaction(const toml::Table& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_global_pair_interaction(), 0);
    using real_type = typename traitsT::real_type;

    const auto potential = get_toml_value<std::string>(
            global, "potential", "[[forcefield.global]]");

    if(potential == "ExcludedVolume")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Excluded Volume.");
        using potential_t = ExcludedVolumePotential<real_type>;

        return read_spatial_partition<traitsT, potential_t>(
            global, read_excluded_volume_potential<real_type>(global));
    }
    else if(potential == "DebyeHuckel")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Debye-Huckel.");
        using potential_t = DebyeHuckelPotential<real_type>;

        return read_spatial_partition<traitsT, potential_t>(
            global, read_debye_huckel_potential<real_type>(global));
    }
    else if(potential == "LennardJones")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Lennard-Jones.");
        using potential_t = LennardJonesPotential<real_type>;

        return read_spatial_partition<traitsT, potential_t>(
            global, read_lennard_jones_potential<real_type>(global));
    }
    else if(potential == "UniformLennardJones")
    {
        MJOLNIR_LOG_NOTICE("-- potential functions is Uniform Lennard-Jones.");
        using potential_t = UniformLennardJonesPotential<real_type>;

        return read_spatial_partition<traitsT, potential_t>(
            global, read_uniform_lennard_jones_potential<real_type>(global));
    }
    else
    {
        throw_exception<std::runtime_error>(
                "invalid potential as GlobalPairInteraction: ", potential);
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
#endif// MJOLNIR_READ_INTERACTION
