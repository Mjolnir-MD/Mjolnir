#ifndef MJOLNIR_INPUT_READ_GLOBAL_POTENTIAL_HPP
#define MJOLNIR_INPUT_READ_GLOBAL_POTENTIAL_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/input/utility.hpp>
#include <mjolnir/forcefield/global/ExcludedVolumePotential.hpp>
#include <mjolnir/forcefield/global/InversePowerPotential.hpp>
#include <mjolnir/forcefield/global/HardCoreExcludedVolumePotential.hpp>
#include <mjolnir/forcefield/global/LennardJonesPotential.hpp>
#include <mjolnir/forcefield/global/UniformLennardJonesPotential.hpp>
#include <mjolnir/forcefield/global/DebyeHuckelPotential.hpp>
#include <mjolnir/forcefield/3SPN2/ThreeSPN2ExcludedVolumePotential.hpp>
#include <mjolnir/core/Topology.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// ============================================================================
// global potential
// ============================================================================

//
// reads `ignore.molecule`. If not provided, return the default value.
//
inline IgnoreMolecule<typename Topology::molecule_id_type>
read_ignored_molecule(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    if(global.as_table().count("ignore") == 0)
    {
        MJOLNIR_LOG_NOTICE("No `ignore.molecule` is provided. "
                           "All the molecules are taken into account.");

        return IgnoreMolecule<typename Topology::molecule_id_type>{
            make_unique<IgnoreNothing<typename Topology::molecule_id_type>>()
        };
    }

    const auto& ignore = toml::find(global, "ignore");
    check_keys_available(ignore, {"molecule", "group", "particles_within"});

    if(ignore.as_table().count("molecule") == 0)
    {
        MJOLNIR_LOG_NOTICE("No `ignore.molecule` is provided. "
                           "All the groups are taken into account.");

        return IgnoreMolecule<typename Topology::molecule_id_type>{
            make_unique<IgnoreNothing<typename Topology::molecule_id_type>>()
        };
    }

    const auto name = toml::find<std::string>(ignore, "molecule");

    if(name == "Nothing")
    {
        MJOLNIR_LOG_NOTICE("all the interactions"
                           "(both (inter|intra)-molecule) are included");
        return IgnoreMolecule<typename Topology::molecule_id_type>{
            make_unique<IgnoreNothing<typename Topology::molecule_id_type>>()
        };
    }
    else if(name == "Self" || name == "Intra")
    {
        MJOLNIR_LOG_NOTICE("intra-molecule interaction is ignored");
        return IgnoreMolecule<typename Topology::molecule_id_type>{
            make_unique<IgnoreSelf<typename Topology::molecule_id_type>>()
        };
    }
    else if(name == "Others" || name == "Inter")
    {
        MJOLNIR_LOG_NOTICE("inter-molecule interaction is ignored");
        return IgnoreMolecule<typename Topology::molecule_id_type>{
            make_unique<IgnoreOthers<typename Topology::molecule_id_type>>()
        };
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error(
            "[error] mjolnir::read_ignored_molecule: unknown setting",
            toml::find(ignore, "molecule"), "expected (Nothing|Self|Others)."));
    }
}

//
// reads `ignore.group`. If not provided, return the default value.
//
inline IgnoreGroup<typename Topology::group_id_type>
read_ignored_group(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using group_id_type = typename Topology::group_id_type;

    // ```toml
    // # array of strings
    // ignore.group.intra = ["DNA"] # ignore intra-DNA
    // # pair of strings
    // ignore.group.inter = [
    //     ["DNA", "Protein1"], # ignore DNA-protein1 interaction
    //     ["DNA", "Protein2"]  # ignore DNA-protein2 interaction
    // ]
    // ```

    std::map<group_id_type, std::vector<group_id_type>> ignores;

    if(global.as_table().count("ignore") == 0)
    {
        MJOLNIR_LOG_NOTICE("No `ignore.group` is provided. "
                           "All the groups are taken into account.");
        assert(ignores.empty());
        return IgnoreGroup<group_id_type>(ignores);
    }

    const auto& ignore = toml::find(global, "ignore");

    if(ignore.as_table().count("group") == 0)
    {
        MJOLNIR_LOG_NOTICE("No `ignore.group` is provided. "
                           "All the groups are taken into account.");
        assert(ignores.empty());
        return IgnoreGroup<group_id_type>(ignores);
    }

    const auto group = toml::find(ignore, "group");
    check_keys_available(group, {"intra", "inter"});

    if(group.as_table().count("intra") == 1)
    {
        for(auto intra : toml::find<std::vector<std::string>>(group, "intra"))
        {
            assert(ignores.count(intra) == 0);
            ignores[intra] = {intra};
            MJOLNIR_LOG_NOTICE("ignore interactions inside ", intra);
        }
    }
    if(group.as_table().count("inter") == 1)
    {
        for(auto inter : toml::find<
            std::vector<std::pair<std::string, std::string>>>(group, "inter"))
        {
            const auto fst = std::move(inter.first);
            const auto snd = std::move(inter.second);
            if(ignores.count(fst) == 0) {ignores[fst] = {};}
            if(ignores.count(snd) == 0) {ignores[snd] = {};}

            ignores.at(fst).push_back(snd);
            ignores.at(snd).push_back(fst);

            MJOLNIR_LOG_NOTICE("ignore interactions between ", fst, " and ", snd);
        }
    }
    return IgnoreGroup<group_id_type>(ignores);
}

//
// reads `ignore.particles_within`. If not provided, return the default value.
//
inline std::map<std::string, std::size_t>
read_ignore_particles_within(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();

    using map_type = std::map<std::string, std::size_t>;

    if(global.as_table().count("ignore") == 0)
    {
        return map_type{};
    }
    const auto& ignore = toml::find(global, "ignore");

    const auto ignore_particle_within = toml::find_or<map_type>(
            ignore, "particles_within", map_type{});
    for(const auto& connection : ignore_particle_within)
    {
        MJOLNIR_LOG_NOTICE("particles that have connection ", connection.first,
                           " within ", connection.second, " will be ignored");
    }
    return ignore_particle_within;
}

//
// checks index overlap in the `parameters` array. If any, show the warning.
// After reading the parameters, call this.
// ```toml
// parameters = [
//     {index = 10, radius = 2.0},
//     {index = 10, radius = 2.0}, # <- overlap! parameter value is ambiguous.
// ]
// ```
template<typename parameterT>
void check_parameter_overlap(const toml::value& env, const toml::array& setting,
        std::vector<std::pair<std::size_t, parameterT>>& parameters)
{
    if(parameters.empty()) {return ;}

    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using value_type = std::pair<std::size_t, parameterT>;

    std::sort(parameters.begin(), parameters.end(),
            [](const value_type& lhs, const value_type& rhs) noexcept -> bool {
                return lhs.first < rhs.first; // check only its index.
            });

    const auto overlap = std::adjacent_find(parameters.begin(), parameters.end(),
            [](const value_type& lhs, const value_type& rhs) noexcept -> bool {
                return lhs.first == rhs.first;
            });

    if(overlap != parameters.end())
    {
        const std::size_t overlapped_idx = overlap->first;
        MJOLNIR_LOG_ERROR("parameter for ", overlapped_idx, " defined twice");

        // define a functor to find the toml::value by its index.
        // read the index in parameter and compare it with the `overlapped_idx`.
        const auto find_by_index =
            [overlapped_idx, &env](const toml::value& v) -> bool {
                return find_parameter<std::size_t>(v, env, "index") +
                    find_parameter_or<std::size_t>(v, env, "offset", 0) ==
                    overlapped_idx;
            };

        const auto overlapped1 = std::find_if(
                setting.begin(), setting.end(), find_by_index);
        assert(overlapped1 != setting.end());

        const auto overlapped2 = std::find_if(
                std::next(overlapped1), setting.end(), find_by_index);
        assert(overlapped2 != setting.end());

        throw_exception<std::runtime_error>(toml::format_error(
            "check_parameter_overlap: duplicate parameter definitions",
            *overlapped1, "this defined twice", *overlapped2, "here"));
    }
    return ;
}

template<typename traitsT>
ExcludedVolumePotential<traitsT>
read_excluded_volume_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using potential_type = ExcludedVolumePotential<traitsT>;
    using real_type      = typename potential_type::real_type;
    using parameter_type = typename potential_type::parameter_type;

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    const real_type eps = toml::find<real_type>(global, "epsilon");
    MJOLNIR_LOG_INFO("epsilon = ", eps);

    const real_type cutoff = toml::find_or<real_type>(global, "cutoff",
            potential_type::default_cutoff());
    MJOLNIR_LOG_INFO("relative cutoff = ", cutoff);

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx    = find_parameter<std::size_t>(param, env, "index") +
                            find_parameter_or<std::size_t>(param, env, "offset", 0);
        const auto radius = find_parameter<real_type  >(param, env, "radius");

        params.emplace_back(idx, radius);
        MJOLNIR_LOG_INFO("idx = ", idx, ", radius = ", radius);
    }
    check_parameter_overlap(env, ps, params);

    return potential_type(eps, cutoff, params,
            read_ignore_particles_within(global),
            read_ignored_molecule(global), read_ignored_group(global));
}

template<typename traitsT>
InversePowerPotential<traitsT>
read_inverse_power_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using potential_type = InversePowerPotential<traitsT>;
    using real_type      = typename potential_type::real_type;
    using integer_type   = typename potential_type::integer_type;
    using parameter_type = typename potential_type::parameter_type;

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    const real_type eps = toml::find<real_type>(global, "epsilon");
    MJOLNIR_LOG_INFO("epsilon = ", eps);

    const integer_type n = toml::find<integer_type>(global, "n");
    MJOLNIR_LOG_INFO("n = ", n);

    const real_type cutoff = toml::find_or<real_type>(global, "cutoff",
            potential_type::default_cutoff(n));
    MJOLNIR_LOG_INFO("relative cutoff = ", cutoff);

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx    = find_parameter<std::size_t>(param, env, "index") +
                            find_parameter_or<std::size_t>(param, env, "offset", 0);
        const auto radius = find_parameter<real_type  >(param, env, "radius");

        params.emplace_back(idx, radius);
        MJOLNIR_LOG_INFO("idx = ", idx, ", radius = ", radius);
    }
    check_parameter_overlap(env, ps, params);

    return potential_type(eps, n, cutoff, params,
            read_ignore_particles_within(global),
            read_ignored_molecule(global), read_ignored_group(global));
}

template<typename traitsT>
HardCoreExcludedVolumePotential<traitsT>
read_hard_core_excluded_volume_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using potential_type = HardCoreExcludedVolumePotential<traitsT>;
    using real_type      = typename potential_type::real_type;
    using parameter_type = typename potential_type::parameter_type;

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    const real_type eps = toml::find<real_type>(global, "epsilon");
    MJOLNIR_LOG_INFO("epsilon = ", eps);

    const real_type cutoff = toml::find_or<real_type>(global, "cutoff",
            potential_type::default_cutoff());
    MJOLNIR_LOG_INFO("relative cutoff = ", cutoff);

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx = find_parameter<std::size_t>(param, env, "index") +
                         find_parameter_or<std::size_t>(param, env, "offset", 0);

        const auto core_radius          =
            find_parameter<real_type>(param, env, "core_radius");
        const auto soft_shell_thickness =
            find_parameter<real_type>(param, env, "soft_shell_thickness");

        params.emplace_back(idx, parameter_type{soft_shell_thickness, core_radius});
        MJOLNIR_LOG_INFO("idx = ", idx, ", core_radius = ", core_radius,
                         ", soft_shell_thickness = ", soft_shell_thickness);
    }

    check_parameter_overlap(env, ps, params);

    return potential_type(eps, cutoff, std::move(params),
        read_ignore_particles_within(global),
        read_ignored_molecule(global), read_ignored_group(global));
}

template<typename traitsT>
LennardJonesPotential<traitsT>
read_lennard_jones_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using potential_type = LennardJonesPotential<traitsT>;
    using real_type      = typename potential_type::real_type;
    using parameter_type = typename potential_type::parameter_type;

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    const real_type cutoff = toml::find_or<real_type>(global, "cutoff",
            potential_type::default_cutoff());
    MJOLNIR_LOG_INFO("relative cutoff = ", cutoff);

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx     = find_parameter<std::size_t>(param, env, "index") +
                             find_parameter_or<std::size_t>(param, env, "offset", 0);
        const auto sigma   = find_parameter<real_type>(param, env, "sigma",   u8"σ");
        const auto epsilon = find_parameter<real_type>(param, env, "epsilon", u8"ε");

        params.emplace_back(idx, parameter_type{sigma, epsilon});
        MJOLNIR_LOG_INFO("idx = ", idx, ", sigma = ", sigma, ", epsilon = ", epsilon);
    }

    check_parameter_overlap(env, ps, params);

    return potential_type(cutoff, std::move(params),
            read_ignore_particles_within(global),
            read_ignored_molecule(global), read_ignored_group(global));
}

template<typename traitsT>
UniformLennardJonesPotential<traitsT>
read_uniform_lennard_jones_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using potential_type = UniformLennardJonesPotential<traitsT>;
    using real_type      = typename potential_type::real_type;
    using parameter_type = typename potential_type::parameter_type;

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    const auto sigma   = toml::expect<real_type>(global, u8"σ").or_other(
                         toml::expect<real_type>(global, "sigma")).unwrap();
    const auto epsilon = toml::expect<real_type>(global, u8"ε").or_other(
                         toml::expect<real_type>(global, "epsilon")).unwrap();
    MJOLNIR_LOG_INFO("sigma   = ", sigma);
    MJOLNIR_LOG_INFO("epsilon = ", epsilon);

    const real_type cutoff = toml::find_or<real_type>(global, "cutoff",
            potential_type::default_cutoff());
    MJOLNIR_LOG_INFO("relative cutoff = ", cutoff);

    std::vector<std::pair<std::size_t, parameter_type>> params;
    if(global.as_table().count("parameters") == 1)
    {
        const auto& parameters = toml::find<toml::array>(global, "parameters");
        for(const auto& param : parameters)
        {
            const auto idx = find_parameter<std::size_t>(param, env, "index") +
                             find_parameter_or<std::size_t>(param, env, "offset", 0);
            params.emplace_back(idx, parameter_type{});
        }
        check_parameter_overlap(env, parameters, params);
    }
    else
    {
        MJOLNIR_LOG_WARN("deprecated: `parameters` field in UniformLennardJones"
                         "potential will be required in the future release.");
        MJOLNIR_LOG_WARN("deprecated: write participants explicitly.");
    }

    return potential_type(sigma, epsilon, cutoff, params,
            read_ignore_particles_within(global),
            read_ignored_molecule(global), read_ignored_group(global));
}

template<typename traitsT>
DebyeHuckelPotential<traitsT>
read_debye_huckel_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using potential_type = DebyeHuckelPotential<traitsT>;
    using real_type      = typename potential_type::real_type;
    using parameter_type = typename potential_type::parameter_type;

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    const real_type cutoff = toml::find_or<real_type>(global, "cutoff",
            potential_type::default_cutoff());
    MJOLNIR_LOG_INFO("relative cutoff = ", cutoff);

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx = find_parameter<std::size_t>(param, env, "index") +
                         find_parameter_or<std::size_t>(param, env, "offset", 0);
        const auto charge = find_parameter<real_type  >(param, env, "charge");

        params.emplace_back(idx, parameter_type{charge});
        MJOLNIR_LOG_INFO("idx = ", idx, ", charge = ", charge);
    }

    check_parameter_overlap(env, ps, params);

    return potential_type(cutoff, std::move(params),
            read_ignore_particles_within(global),
            read_ignored_molecule(global), read_ignored_group(global));
}

template<typename traitsT>
ThreeSPN2ExcludedVolumePotential<traitsT>
read_3spn2_excluded_volume_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using potential_type = ThreeSPN2ExcludedVolumePotential<traitsT>;
    using real_type      = typename potential_type::real_type;
    using parameter_type = typename potential_type::parameter_type;
    using bead_kind      = parameter_3SPN2::bead_kind;

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx = find_parameter<std::size_t>(param, env, "index") +
                         find_parameter_or<std::size_t>(param, env, "offset", 0);

        const auto kind = toml::find<std::string>(param, "kind");
        if(kind != "S" && kind != "P" &&
           kind != "A" && kind != "T" && kind != "G" && kind != "C")
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_3spn2_excluded_volume_potential: unknown bead "
                "kind", find_parameter<toml::value>(param, env, "kind"),
                "expected S, P, A, T, G, C."));
        }
        bead_kind bead;
        switch(kind.front())
        {
            case 'P': {bead = bead_kind::Phosphate; break;}
            case 'S': {bead = bead_kind::Sugar;     break;}
            case 'A': {bead = bead_kind::BaseA;     break;}
            case 'T': {bead = bead_kind::BaseT;     break;}
            case 'G': {bead = bead_kind::BaseG;     break;}
            case 'C': {bead = bead_kind::BaseC;     break;}
            default:  {assert(false);}
        }
        params.emplace_back(idx, bead);
        MJOLNIR_LOG_INFO("idx = ", idx, ", kind = ", bead);
    }

    ThreeSPN2ExcludedVolumePotentialParameter<real_type> default_parameters;

    check_parameter_overlap(env, ps, params);
    return potential_type(std::move(default_parameters), params,
            read_ignore_particles_within(global),
            read_ignored_molecule(global), read_ignored_group(global));
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template ExcludedVolumePotential<SimulatorTraits<double, UnlimitedBoundary>       > read_excluded_volume_potential(const toml::value& global);
extern template ExcludedVolumePotential<SimulatorTraits<float,  UnlimitedBoundary>       > read_excluded_volume_potential(const toml::value& global);
extern template ExcludedVolumePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_excluded_volume_potential(const toml::value& global);
extern template ExcludedVolumePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_excluded_volume_potential(const toml::value& global);

extern template InversePowerPotential<SimulatorTraits<double, UnlimitedBoundary>       > read_inverse_power_potential(const toml::value&);
extern template InversePowerPotential<SimulatorTraits<float,  UnlimitedBoundary>       > read_inverse_power_potential(const toml::value&);
extern template InversePowerPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_inverse_power_potential(const toml::value&);
extern template InversePowerPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_inverse_power_potential(const toml::value&);

extern template HardCoreExcludedVolumePotential<SimulatorTraits<double, UnlimitedBoundary>       > read_hard_core_excluded_volume_potential(const toml::value& global);
extern template HardCoreExcludedVolumePotential<SimulatorTraits<float,  UnlimitedBoundary>       > read_hard_core_excluded_volume_potential(const toml::value& global);
extern template HardCoreExcludedVolumePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_hard_core_excluded_volume_potential(const toml::value& global);
extern template HardCoreExcludedVolumePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_hard_core_excluded_volume_potential(const toml::value& global);

extern template LennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>       > read_lennard_jones_potential(const toml::value& global);
extern template LennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>       > read_lennard_jones_potential(const toml::value& global);
extern template LennardJonesPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_lennard_jones_potential(const toml::value& global);
extern template LennardJonesPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_lennard_jones_potential(const toml::value& global);

extern template UniformLennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>       > read_uniform_lennard_jones_potential(const toml::value& global);
extern template UniformLennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>       > read_uniform_lennard_jones_potential(const toml::value& global);
extern template UniformLennardJonesPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_uniform_lennard_jones_potential(const toml::value& global);
extern template UniformLennardJonesPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_uniform_lennard_jones_potential(const toml::value& global);

extern template DebyeHuckelPotential<SimulatorTraits<double, UnlimitedBoundary>       > read_debye_huckel_potential(const toml::value& global);
extern template DebyeHuckelPotential<SimulatorTraits<float,  UnlimitedBoundary>       > read_debye_huckel_potential(const toml::value& global);
extern template DebyeHuckelPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_debye_huckel_potential(const toml::value& global);
extern template DebyeHuckelPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_debye_huckel_potential(const toml::value& global);

extern template ThreeSPN2ExcludedVolumePotential<SimulatorTraits<double, UnlimitedBoundary>       > read_3spn2_excluded_volume_potential(const toml::value& global);
extern template ThreeSPN2ExcludedVolumePotential<SimulatorTraits<float,  UnlimitedBoundary>       > read_3spn2_excluded_volume_potential(const toml::value& global);
extern template ThreeSPN2ExcludedVolumePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_3spn2_excluded_volume_potential(const toml::value& global);
extern template ThreeSPN2ExcludedVolumePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_3spn2_excluded_volume_potential(const toml::value& global);
#endif

} // mjolnir
#endif // MJOLNIR_READ_GLOBAL_POTENTIAL_HPP
