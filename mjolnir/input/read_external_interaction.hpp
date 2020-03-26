#ifndef MJOLNIR_INPUT_READ_EXTERNAL_INTERACTION_HPP
#define MJOLNIR_INPUT_READ_EXTERNAL_INTERACTION_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/forcefield/AFMFit/AFMFitInteraction.hpp>
#include <mjolnir/forcefield/external/ExternalDistanceInteraction.hpp>
#include <mjolnir/forcefield/external/PositionRestraintInteraction.hpp>
#include <mjolnir/core/AxisAlignedPlane.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_external_potential.hpp>
#include <mjolnir/input/read_local_potential.hpp>
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
    MJOLNIR_LOG_FUNCTION();
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
            toml::find(external, "potential"), "here", {
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
    MJOLNIR_LOG_FUNCTION();
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
            using shape_t = AxisAlignedPlane<traitsT, PositiveXDirection<traitsT>>;

            return read_external_distance_interaction<traitsT, shape_t>(
                    external, shape_t(pos, margin));
        }
        else if(axis == "-X")
        {
            MJOLNIR_LOG_NOTICE("-- interaction shape is YZ-Plane (-).");
            using shape_t = AxisAlignedPlane<traitsT, NegativeXDirection<traitsT>>;
            return read_external_distance_interaction<traitsT, shape_t>(
                    external, shape_t(pos, margin));
        }
        else if(axis == "Y" || axis == "+Y")
        {
            MJOLNIR_LOG_NOTICE("-- interaction shape is ZX-Plane (+).");
            using shape_t = AxisAlignedPlane<traitsT, PositiveYDirection<traitsT>>;
            return read_external_distance_interaction<traitsT, shape_t>(
                    external, shape_t(pos, margin));
        }
        else if(axis == "-Y")
        {
            MJOLNIR_LOG_NOTICE("-- interaction shape is ZX-Plane (-).");
            using shape_t = AxisAlignedPlane<traitsT, NegativeYDirection<traitsT>>;
            return read_external_distance_interaction<traitsT, shape_t>(
                    external, shape_t(pos, margin));
        }
        else if(axis == "Z" || axis == "+Z")
        {
            MJOLNIR_LOG_NOTICE("-- interaction shape is XY-Plane (+).");
            using shape_t = AxisAlignedPlane<traitsT, PositiveZDirection<traitsT>>;
            return read_external_distance_interaction<traitsT, shape_t>(
                    external, shape_t(pos, margin));
        }
        else if(axis == "-Z")
        {
            MJOLNIR_LOG_NOTICE("-- interaction shape is XY-Plane (-).");
            using shape_t = AxisAlignedPlane<traitsT, NegativeZDirection<traitsT>>;
            return read_external_distance_interaction<traitsT, shape_t>(
                    external, shape_t(pos, margin));
        }
        else
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_external_distance_interaction_shape: invalid shape",
                toml::find(external, "shape"), "expected [+-]?[XYZ]"));
        }
    }
    else
    {
        throw std::runtime_error("invalid external forcefield shape: " + name);
    }
}

template<typename traitsT>
std::unique_ptr<ExternalForceInteractionBase<traitsT>>
read_external_position_restraint_interaction(const toml::value& external)
{
    // [[forcefields.external]]
    // interaction = "PositionRestraint"
    // potential   = "Harmonic" # ----------------------------------------+
    // parameters  = [                                                    |
    //     {index = 0, position = [0.0, 0.0, 0.0], k = 0.1, v0 = 10.0},   v
    //     #                                       ^^^^^^^^^^^^^^^^^^ Harmonic
    //     # ...
    // ]
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type       = typename traitsT::real_type;
    using coordinate_type = typename traitsT::coordinate_type;

    const auto& env = external.as_table().count("env") == 1 ?
                      external.as_table().at("env") : toml::value{};

    const auto potential = toml::find<std::string>(external, "potential");
    if(potential == "Harmonic")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Harmonic.");
        using potential_t   = HarmonicPotential<real_type>;
        using interaction_t = PositionRestraintInteraction<traitsT, potential_t>;

        // container to put parameters
        typename interaction_t::potential_container_type potentials;

        const auto& parameters = toml::find<toml::array>(external, "parameters");
        potentials.reserve(parameters.size());

        for(auto para : parameters)
        {
            const auto idx = find_parameter<std::size_t>(para, env, "index") +
                             find_parameter_or<std::int64_t>(para, env, "offset", 0);
            const auto pos = find_parameter<std::array<real_type, 3>>(para, env, "position");
            const auto crd = math::make_coordinate<coordinate_type>(pos[0], pos[1], pos[2]);

            // suppress warnings in read_harmonic_potential
            para.as_table().erase("index");
            para.as_table().erase("position");

            const auto pot = read_harmonic_potential<real_type>(para, env);
            potentials.emplace_back(idx, crd, pot);
        }
        return make_unique<interaction_t>(std::move(potentials));
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_external_position_restraint_interaction: invalid potential",
            toml::find(external, "potential"), "here", {
            "expected value is one of the following.",
            "- \"Harmonic\": very famous `k(x-x0)^2` potential"
            }));
    }
}

template<typename traitsT>
std::unique_ptr<ExternalForceInteractionBase<traitsT>>
read_afm_flexible_fitting_interaction(const toml::value& external)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type     = typename traitsT::real_type;
    using interaction_t = AFMFitInteraction<traitsT>;

    const auto k        = toml::find<real_type  >(external, "k");
    const auto gamma    = toml::find<real_type  >(external, "gamma");
    const auto z0       = toml::find<real_type  >(external, "z0");
    const auto cutoff   = toml::find<real_type  >(external, "cutoff");
    const auto margin   = toml::find<real_type  >(external, "margin");
    const auto sigma_x  = toml::find<real_type  >(external, "sigma_x");
    const auto sigma_y  = toml::find<real_type  >(external, "sigma_y");
    const auto pixel_x  = toml::find<real_type  >(external, "pixel_x");
    const auto pixel_y  = toml::find<real_type  >(external, "pixel_y");
    const auto length_x = toml::find<std::size_t>(external, "length_x");
    const auto length_y = toml::find<std::size_t>(external, "length_y");

    const auto& env = external.as_table().count("env") != 0 ?
                      external.at("env") : toml::value{};

    // read radius and participants

    const auto& parameters = toml::find<toml::array>(external, "parameters");

    std::vector<std::pair<std::size_t, real_type>> radii;
    radii.reserve(parameters.size());
    for(const auto& para : parameters)
    {
        const auto idx = find_parameter<std::size_t>(para, env, "index");
        const auto rad = find_parameter<real_type  >(para, env, "radius");
        radii.emplace_back(idx, rad);
    }

    std::vector<real_type> image;
    image.reserve(length_x * length_y);
    for(real_type img : toml::find<std::vector<real_type>>(external, "image"))
    {
        image.push_back(img);
    }
    return make_unique<interaction_t>(k, gamma, z0, cutoff, margin,
            sigma_x, sigma_y, pixel_x, pixel_y, length_x, length_y,
            std::move(radii), std::move(image));
}

// ----------------------------------------------------------------------------
// general read_external_interaction function
// ----------------------------------------------------------------------------

template<typename traitsT>
std::unique_ptr<ExternalForceInteractionBase<traitsT>>
read_external_interaction(const toml::value& external)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    const auto interaction = toml::find<std::string>(external, "interaction");

    if(interaction == "Distance")
    {
        MJOLNIR_LOG_NOTICE("Distance interaction found.");
        return read_external_distance_interaction_shape<traitsT>(external);
    }
    else if(interaction == "PositionRestraint")
    {
        MJOLNIR_LOG_NOTICE("PositionRestraint interaction found.");
        return read_external_position_restraint_interaction<traitsT>(external);
    }
    else if(interaction == "AFMFlexibleFitting")
    {
        MJOLNIR_LOG_NOTICE("AFMFlexibleFitting interaction found.");
        return read_afm_flexible_fitting_interaction<traitsT>(external);
    }
    else
    {
        throw std::runtime_error(toml::format_error("[error] "
            "mjolnir::read_external_interaction: invalid interaction",
            toml::find(external, "interaction"), "here", {
            "expected value is one of the following.",
            "- \"Distance\":           interaction depending on the distance between a particle and a structure",
            "- \"PositionRestraint\":  interaction depending on the distance between a particle and a fixed point",
            "- \"AFMFlexibleFitting\": interaction that fits molecule to an AFM image"
            }));
    }
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_external_interaction(const toml::value& external);
extern template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_external_interaction(const toml::value& external);
extern template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_external_interaction(const toml::value& external);
extern template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_external_interaction(const toml::value& external);

extern template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_external_position_restraint_interaction(const toml::value& external);
extern template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_external_position_restraint_interaction(const toml::value& external);
extern template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_external_position_restraint_interaction(const toml::value& external);
extern template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_external_position_restraint_interaction(const toml::value& external);

extern template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_external_distance_interaction_shape(const toml::value& external);
extern template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_external_distance_interaction_shape(const toml::value& external);
extern template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_external_distance_interaction_shape(const toml::value& external);
extern template std::unique_ptr<ExternalForceInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_external_distance_interaction_shape(const toml::value& external);
#endif

} // mjolnir
#endif// MJOLNIR_READ_EXTERNAL_INTERACTION_HPP
