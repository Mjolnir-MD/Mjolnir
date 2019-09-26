#ifndef MJOLNIR_INPUT_READ_LOCAL_INTERACTION_HPP
#define MJOLNIR_INPUT_READ_LOCAL_INTERACTION_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/interaction/local/BondLengthInteraction.hpp>
#include <mjolnir/interaction/local/ContactInteraction.hpp>
#include <mjolnir/interaction/local/DirectionalContactInteraction.hpp>
#include <mjolnir/interaction/local/BondAngleInteraction.hpp>
#include <mjolnir/interaction/local/DihedralAngleInteraction.hpp>
#include <mjolnir/interaction/local/DummyInteraction.hpp>

#include <mjolnir/forcefield/3SPN2/ThreeSPN2BaseStackingPotential.hpp>
#include <mjolnir/forcefield/3SPN2/ThreeSPN2BaseStackingInteraction.hpp>

#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_local_potential.hpp>
#include <memory>
#include <utility>

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
    else if(potential == "3SPN2Bond")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is 3SPN2Bond.");
        using potentialT = ThreeSPN2BondPotential<real_type>;

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
            "- \"Gaussian\" : well-known gaussian potential",
            "- \"3SPN2Bond\": bond length potential for 3SPN2"
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

template<std::size_t N, typename realT, typename angle1_potentialT,
         typename angle2_potentialT, typename contact_potentialT>
std::vector<std::tuple<std::array<std::size_t, N>, angle1_potentialT,
                       angle2_potentialT, contact_potentialT>>
read_directional_contact_potentials(const toml::value& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    MJOLNIR_LOG_INFO("as", N, "-body interaction");

    using indices_type = std::array<std::size_t, N>;
    using indices_potentials_tuple_type = std::tuple<
        indices_type, angle1_potentialT, angle2_potentialT, contact_potentialT
        >;

    const auto& params = toml::find<toml::array>(local, "parameters");
    MJOLNIR_LOG_NOTICE("-- ", params.size(), " interactions are found.");

    const auto& env = local.as_table().count("env") == 1 ?
            local.as_table().at("env") : toml::value{};

    std::vector<indices_potentials_tuple_type> retval;
    retval.reserve(params.size());
    for(const auto& item : params)
    {
        const auto indices = find_parameter<indices_type>(item, env, "indices");
        MJOLNIR_LOG_INFO_NO_LF("idxs = ", indices, ", ");

        const auto angle_pot1 = find_parameter<toml::value>(item, env, "angle1");
        const auto angle_pot2 = find_parameter<toml::value>(item, env, "angle2");
        const auto contact_pot = find_parameter<toml::value>(item, env, "contact");

        retval.emplace_back(
            std::make_tuple(indices,
                detail::read_local_potential_impl<angle1_potentialT>::invoke(angle_pot1, env),
                detail::read_local_potential_impl<angle2_potentialT>::invoke(angle_pot2, env),
                detail::read_local_potential_impl<contact_potentialT>::invoke(contact_pot, env))
            );
    }
    return retval;
}

// This is reading contact part of read_directional_contact_interaction function.
template<typename traitsT, typename ... PotentialTs>
typename std::enable_if<
    sizeof...(PotentialTs) == 2, std::unique_ptr<LocalInteractionBase<traitsT>>
    >::type
read_directional_contact_interaction(
    const std::string& kind, const toml::value& local, std::vector<std::string>&&)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = typename traitsT::real_type;

    real_type margin = 0.5; // default value
    if (local.as_table().count("margin") == 1)
    {
        margin = toml::find<real_type>(local, "margin");
    }

    const auto contact_potential = toml::find<std::string>(local, "potentials", "contact");
    if(contact_potential == "GoContact")
    {
        MJOLNIR_LOG_NOTICE("-- contact potential function is 10-12 Go contact.");
        using contact_potentialT = GoContactPotential<real_type>;

        return make_unique<
            DirectionalContactInteraction<
                traitsT, PotentialTs..., contact_potentialT
                >
            >(kind,
              read_directional_contact_potentials<
                  4, real_type, PotentialTs..., contact_potentialT
              >(local), margin);
    }
    else if(contact_potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- contact potential function is Gaussian.");
        using contact_potentialT = GaussianPotential<real_type>;

        return make_unique<
            DirectionalContactInteraction<
                traitsT, PotentialTs..., contact_potentialT
                >
            >(kind,
              read_directional_contact_potentials<
                  4, real_type, PotentialTs..., contact_potentialT
              >(local), margin);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_directional_contact_interaction: invalid contact potential",
            toml::find<toml::value>(local, "potentials", "contact"), "here", {
            "expected value is one of the following.",
            "- \"GoContact\": r^12 - r^10 type native contact potential",
            "- \"Gaussian\" : well-known gaussian potential"
            }));
    }
}

// This is reading angle parts of read_directional_contact_interaction function.
template<typename traitsT, typename ... PotentialTs>
typename std::enable_if<
    sizeof...(PotentialTs) < 2, std::unique_ptr<LocalInteractionBase<traitsT>>
    >::type
read_directional_contact_interaction(
    const std::string& kind, const toml::value& local, std::vector<std::string>&& angle_potential_keys)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = typename traitsT::real_type;

    const std::string angle_potential_key = angle_potential_keys.back();
    const auto        angle_potential = toml::find<std::string>(local, "potentials", angle_potential_key);
    angle_potential_keys.pop_back();
    if(angle_potential == "Cosine")
    {
        MJOLNIR_LOG_NOTICE("-- angle potential function is Cosine potential");
        using angle_potential_T = CosinePotential<real_type>;

        return read_directional_contact_interaction<
            traitsT, PotentialTs..., angle_potential_T
            >(kind, local, std::move(angle_potential_keys));
    }
    else if(angle_potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- angle potential function is Gaussian");
        using angle_potential_T = GaussianPotential<real_type>;

        return read_directional_contact_interaction<
            traitsT, PotentialTs..., angle_potential_T
            >(kind, local, std::move(angle_potential_keys));
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_bond_length_interaction: invalid angle potential",
            toml::find<toml::value>(local, "potentials", angle_potential_key), "here", {
            "expected value is one of the following.",
            "- \"Cosine\"   : 1 + Cosine(x) potential",
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
    using real_type      = typename traitsT::real_type;
    using potential_type = ThreeSPN2BaseStackingPotential<real_type>;
    using base_kind      = parameter_3SPN2::base_kind;

    // [[forcefields.local]]
    // interaction = "3SPN2BaseStacking"
    // topology    = "none"
    // parameters = [
    //     {S_idx = 0, B3_idx = 1, B5_idx = 4, Base3 = "A", Base5 = "T"},
    //     # ...
    // ]

    using indices_type   = std::array<std::size_t, 3>; // {S, B3, B5}
    using parameter_type = typename potential_type::parameter_type;

    potential_type potential;

    const auto& params = toml::find<toml::array>(local, "parameters");
    MJOLNIR_LOG_NOTICE("-- ", params.size(), " interactions are found.");

    const auto& env = local.as_table().count("env") == 1 ?
                      local.as_table().at("env") : toml::value{};

    std::vector<std::pair<indices_type, parameter_type>> parameters;
    parameters.reserve(params.size());
    for(const auto& item : params)
    {
        const auto S_idx  = find_parameter<std::size_t>(item, env, "S_idx");
        const auto B5_idx = find_parameter<std::size_t>(item, env, "B5_idx");
        const auto B3_idx = find_parameter<std::size_t>(item, env, "B3_idx");

        const indices_type idxs{{S_idx, B5_idx, B3_idx}};

        const auto B5 = find_parameter<std::string>(item, env, "Base5");
        const auto B3 = find_parameter<std::string>(item, env, "Base3");
        if(B5 != "A" && B5 != "T" && B5 != "G" && B5 != "C")
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_3spn2_base_stacking_potential: none of A,T,C,G",
                find_parameter<toml::value>(item, env, "Base5"), "here"));
        }
        if(B3 != "A" && B3 != "T" && B3 != "G" && B3 != "C")
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_3spn2_base_stacking_potential: none of A,T,C,G",
                find_parameter<toml::value>(item, env, "Base3"), "here"));
        }
        base_kind B3_kind = base_kind::X;
        base_kind B5_kind = base_kind::X;
        switch(B5.front())
        {
            case 'A': {B5_kind = base_kind::A; break;}
            case 'T': {B5_kind = base_kind::T; break;}
            case 'G': {B5_kind = base_kind::G; break;}
            case 'C': {B5_kind = base_kind::C; break;}
        }
        switch(B3.front())
        {
            case 'A': {B3_kind = base_kind::A; break;}
            case 'T': {B3_kind = base_kind::T; break;}
            case 'G': {B3_kind = base_kind::G; break;}
            case 'C': {B3_kind = base_kind::C; break;}
        }
        assert(B3_kind != base_kind::X);
        assert(B5_kind != base_kind::X);

        const auto bs_kind = potential.bs_kind(B5_kind, B3_kind);

        MJOLNIR_LOG_INFO("ThreeSPN2BaseStackingPotential = {S = ", S_idx,
            ", B5 = ", B5_idx, ", B3 = ", B3_idx, ", bases = ", bs_kind, '}');

        parameters.emplace_back(idxs, bs_kind);
    }
    return make_unique<ThreeSPN2BaseStackingInteraction<traitsT>>(
            kind, std::move(parameters), std::move(potential));
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
    else if(interaction == "DirectionalContact")
    {
        MJOLNIR_LOG_NOTICE("Directional Contact interaction found.");
        return read_directional_contact_interaction<traitsT>(kind, local, {"angle2", "angle1"});
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
            "- \"BondLength\"         : 2-body well-known chemical bond interaction",
            "- \"DirectionalContact\" : 4-body contact interaction depends on the contact angle",
            "- \"BondAngle\"          : 3-body well-known bond angle interaction",
            "- \"DihedralAngle\"      : 4-body well-known dihedral angle interaction",
            "- \"Dummy\"              : To represent a strange topology. It does nothing",
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

extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_3spn2_base_stacking_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_3spn2_base_stacking_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_3spn2_base_stacking_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_3spn2_base_stacking_interaction(const std::string& kind, const toml::value& local);
#endif


} // mjolnir
#endif// MJOLNIR_READ_LOCAL_INTERACTION
