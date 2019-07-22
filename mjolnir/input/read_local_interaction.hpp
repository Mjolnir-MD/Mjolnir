#ifndef MJOLNIR_INPUT_READ_LOCAL_INTERACTION_HPP
#define MJOLNIR_INPUT_READ_LOCAL_INTERACTION_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/interaction/local/BondLengthInteraction.hpp>
#include <mjolnir/interaction/local/ContactInteraction.hpp>
#include <mjolnir/interaction/local/BondAngleInteraction.hpp>
#include <mjolnir/interaction/local/DihedralAngleInteraction.hpp>
#include <mjolnir/interaction/local/DummyInteraction.hpp>
#include <mjolnir/interaction/local/ThreeSPN2BaseStackingInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_local_potential.hpp>
#include <memory>

namespace mjolnir
{

// ----------------------------------------------------------------------------
// individual local interaction. would be called from read_local_interaction
// defined at the bottom of this file.
// ----------------------------------------------------------------------------

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_bond_length_interaction(const std::string& kind, const toml::value& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = typename traitsT::real_type;

    const auto potential = toml::find<std::string>(local, "potential");
    if(potential == "Harmonic")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Harmonic.");
        using potentialT = HarmonicPotential<real_type>;

        return make_unique<BondLengthInteraction<traitsT, potentialT>>(
                kind, read_local_potential<2, potentialT>(local));
    }
    else if(potential == "GoContact")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is 10-12 Go contact.");
        using potentialT = GoContactPotential<real_type>;

        return make_unique<BondLengthInteraction<traitsT, potentialT>>(
                kind, read_local_potential<2, potentialT>(local));
    }
    else if(potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Gaussian.");
        using potentialT = GaussianPotential<real_type>;

        return make_unique<BondLengthInteraction<traitsT, potentialT>>(
                kind, read_local_potential<2, potentialT>(local));
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_bond_length_interaction: invalid potential",
            toml::find<toml::value>(local, "potential"), "here", {
            "expected value is one of the following.",
            "- \"Harmonic\" : well-known harmonic potential",
            "- \"GoContact\": r^12 - r^10 type native contact potential",
            "- \"Gaussian\" : well-known gaussian potential"
            }));
    }
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_contact_interaction(const std::string& kind, const toml::value& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = typename traitsT::real_type;

    real_type margin = 0.5; // default value
    if(local.as_table().count("margin") == 1)
    {
        margin = toml::find<real_type>(local, "margin");
    }

    const auto potential = toml::find<std::string>(local, "potential");
    if(potential == "GoContact")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is 10-12 Go contact.");
        using potentialT = GoContactPotential<real_type>;

        return make_unique<ContactInteraction<traitsT, potentialT>>(
                kind, read_local_potential<2, potentialT>(local), margin);
    }
    else if(potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Gaussian.");
        using potentialT = GaussianPotential<real_type>;

        return make_unique<ContactInteraction<traitsT, potentialT>>(
                kind, read_local_potential<2, potentialT>(local), margin);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_bond_length_interaction: invalid potential",
            toml::find<toml::value>(local, "potential"), "here", {
            "expected value is one of the following.",
            "- \"GoContact\": r^12 - r^10 type native contact potential",
            "- \"Gaussian\" : well-known gaussian potential"
            }));
    }
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_bond_angle_interaction(const std::string& kind, const toml::value& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = typename traitsT::real_type;

    const auto potential = toml::find<std::string>(local, "potential");

    if(potential == "Harmonic")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Harmonic.");
        using potentialT = HarmonicPotential<real_type>;

        return make_unique<BondAngleInteraction<traitsT, potentialT>>(
                kind, read_local_potential<3, potentialT>(local));
    }
    else if(potential == "FlexibleLocalAngle")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Flexible Local Angle.");
        using potentialT = FlexibleLocalAnglePotential<real_type>;

        return make_unique<BondAngleInteraction<traitsT, potentialT>>(
                kind, read_local_potential<3, potentialT>(local));
    }
    else if(potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Gaussian.");
        using potentialT = GaussianPotential<real_type>;

        return make_unique<BondAngleInteraction<traitsT, potentialT>>(
                kind, read_local_potential<3, potentialT>(local));
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_bond_angle_interaction: invalid potential",
            toml::find<toml::value>(local, "potential"), "here", {
            "expected value is one of the following.",
            "- \"Harmonic\"          : well-known harmonic potential",
            "- \"Gaussian\"          : well-known gaussian potential"
            "- \"FlexibleLocalAngle\": table-based potential for C-alpha protein model",
            }));
    }
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_dihedral_angle_interaction(const std::string& kind, const toml::value& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = typename traitsT::real_type;

    const auto potential = toml::find<std::string>(local, "potential");

    if(potential == "Harmonic")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Harmonic.");
        MJOLNIR_LOG_WARN  ("[deprecated] There is a indifferentiable point.");
        MJOLNIR_LOG_WARN  ("[deprecated] Reconsider to use different stuff.");
        using potentialT = HarmonicPotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    else if(potential == "ClementiDihedral")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Clementi-Go's dihedral.");
        using potentialT = ClementiDihedralPotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    else if(potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Gaussian.");
        using potentialT = PeriodicGaussianPotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    else if(potential == "PeriodicGaussian")
    {
        // XXX remove this after v1.4.0
        MJOLNIR_LOG_NOTICE("-- potential function is PeriodicGaussian.");
        MJOLNIR_LOG_WARN  ("[deprecated] \"Periodic\" is no longer needed.");
        using potentialT = PeriodicGaussianPotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    else if(potential == "FlexibleLocalDihedral")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Flexible Local Dihedral.");
        using potentialT = FlexibleLocalDihedralPotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    else if(potential == "Cosine")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Cosine.");
        using potentialT = CosinePotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    else if(potential == "Gaussian+FlexibleLocalDihedral")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is "
                           "Gaussian + FlexibleLocalDihedral.");
        using potentialT = SumLocalPotential<real_type,
              PeriodicGaussianPotential, FlexibleLocalDihedralPotential>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(kind,
            read_local_potentials<4, real_type,
                PeriodicGaussianPotential, FlexibleLocalDihedralPotential>(local));
    }
    else if(potential == "PeriodicGaussian+FlexibleLocalDihedral")
    {
        // XXX remove this after v1.4.0
        MJOLNIR_LOG_NOTICE("-- potential function is "
                           "PeriodicGaussian + FlexibleLocalDihedral.");
        MJOLNIR_LOG_WARN  ("[deprecated] \"Periodic\" is no longer needed.");
        using potentialT = SumLocalPotential<real_type,
              PeriodicGaussianPotential, FlexibleLocalDihedralPotential>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(kind,
            read_local_potentials<4, real_type,
                PeriodicGaussianPotential, FlexibleLocalDihedralPotential>(local));
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_dihedral_angle_interaction: invalid potential",
            toml::find<toml::value>(local, "potential"), "here", {
            "expected value is one of the following.",
            "- \"Harmonic\"             : well-known harmonic potential",
            "- \"Gaussian\"             : well-known gaussian potential"
            "- \"ClementiDihedral\"     : potential used in the off-lattice Go protein model"
            "- \"FlexibleLocalDihedral\": table-based potential for C-alpha protein model",
            }));
    }
}


template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_dummy_interaction(const std::string& kind, const toml::value& local)
{
    // It does not require `potential` field because this interaction does not
    // calculate force or energy.
    //
    // ```toml
    // [[forcefields.local]]
    // interaction = "Dummy"
    // topology  = "bond"
    // parameters = [
    //     {indices = [0, 1]},
    //     # ...
    // ]
    // ```
    //
    // So, unlike other interactions, it reads particle indices here.

    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using indices_type = std::array<std::size_t, 2>;

    if(kind == "none")
    {
        MJOLNIR_LOG_WARN("Dummy potential got `none`. It has no effect.");
        MJOLNIR_LOG_WARN("It is for defining unusual topology.");
    }
    if(local.as_table().count("potential") == 1)
    {
        MJOLNIR_LOG_WARN("Dummy Interaction has a `potential` field.");
        MJOLNIR_LOG_WARN("It is for defining unusual topology.");
        MJOLNIR_LOG_WARN("Potential parameters are ignored.");
    }

    const auto& params = toml::find<toml::array>(local, "parameters");
    MJOLNIR_LOG_NOTICE("-- ", params.size(), " bonds are found.");

    std::vector<indices_type> indices_list;
    indices_list.reserve(params.size());
    for(const auto& item : params)
    {
        const auto indices = toml::find<indices_type>(item, "indices");
        MJOLNIR_LOG_INFO("idxs = ", indices);
        indices_list.push_back(indices);
    }
    return make_unique<DummyInteraction<traitsT>>(kind, std::move(indices_list));
}


template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_3spn2_base_stacking_interaction(const std::string& kind, const toml::value& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type  = typename traitsT::real_type;
    using potentialT = ThreeSPN2BaseStackingPotential<real_type>;

    return make_unique<ThreeSPN2BaseStackingInteraction<traitsT>>(
            kind, read_local_potential<3, potentialT>(local));
}

// ----------------------------------------------------------------------------
// general read_local_interaction function
// ----------------------------------------------------------------------------

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_local_interaction(const toml::value& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto interaction = toml::find<std::string>(local, "interaction");
    const auto kind        = toml::find<std::string>(local, "topology");
    MJOLNIR_LOG_INFO("connection kind = ", kind);

    if(interaction == "BondLength")
    {
        MJOLNIR_LOG_NOTICE("Bond Length interaction found.");
        return read_bond_length_interaction<traitsT>(kind, local);
    }
    else if(interaction == "Contact")
    {
        MJOLNIR_LOG_NOTICE("Contact interaction found.");
        return read_contact_interaction<traitsT>(kind, local);
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
    else if(interaction == "Dummy")
    {
        MJOLNIR_LOG_NOTICE("Dummy interaction found.");
        return read_dummy_interaction<traitsT>(kind, local);
    }
    else if(interaction == "3SPN2BaseStacking")
    {
        MJOLNIR_LOG_NOTICE("3SPN2BaseStacking interaction found.");
        return read_3spn2_base_stacking_interaction<traitsT>(kind, local);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_local_interaction: invalid interaction",
            toml::find<toml::value>(local, "interaction"), "here", {
            "expected value is one of the following.",
            "- \"BondLength\"    : 2-body well-known chemical bond interaction",
            "- \"BondAngle\"     : 3-body well-known bond angle interaction",
            "- \"DihedralAngle\" : 4-body well-known dihedral angle interaction",
            "- \"Dummy\"         : To represent a strange topology. It does nothing",
            }));
    }
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_local_interaction(const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_local_interaction(const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_local_interaction(const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_local_interaction(const toml::value& local);

extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_bond_length_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_bond_length_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_bond_length_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_bond_length_interaction(const std::string& kind, const toml::value& local);

extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_contact_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_contact_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_contact_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_contact_interaction(const std::string& kind, const toml::value& local);

extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_bond_angle_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_bond_angle_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_bond_angle_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_bond_angle_interaction(const std::string& kind, const toml::value& local);

extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_dihedral_angle_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_dihedral_angle_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_dihedral_angle_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_dihedral_angle_interaction(const std::string& kind, const toml::value& local);

extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_dummy_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_dummy_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_dummy_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_dummy_interaction(const std::string& kind, const toml::value& local);
#endif


} // mjolnir
#endif// MJOLNIR_READ_LOCAL_INTERACTION
