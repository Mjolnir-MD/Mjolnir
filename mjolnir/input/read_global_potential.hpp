#ifndef MJOLNIR_INPUT_READ_GLOBAL_POTENTIAL_HPP
#define MJOLNIR_INPUT_READ_GLOBAL_POTENTIAL_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/input/utility.hpp>
#include <mjolnir/forcefield/global/ParameterList.hpp>
#include <mjolnir/forcefield/global/ExcludedVolumePotential.hpp>
#include <mjolnir/forcefield/global/InversePowerPotential.hpp>
#include <mjolnir/forcefield/global/HardCoreExcludedVolumePotential.hpp>
#include <mjolnir/forcefield/global/LennardJonesPotential.hpp>
#include <mjolnir/forcefield/global/LennardJonesAttractivePotential.hpp>
#include <mjolnir/forcefield/global/WCAPotential.hpp>
#include <mjolnir/forcefield/global/DebyeHuckelPotential.hpp>
#include <mjolnir/forcefield/3SPN2/ThreeSPN2ExcludedVolumePotential.hpp>
#include <mjolnir/forcefield/iSoLF/iSoLFAttractivePotential.hpp>
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

template<typename traitsT, typename potentialT, typename Func>
std::unique_ptr<CombinationTable<traitsT, potentialT>>
read_parameter_table(const toml::value& global, Func parameter_reader)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using parameter_list      = CombinationTable<traitsT, potentialT>;
    using parameter_type      = typename parameter_list::parameter_type;
    using pair_parameter_type = typename parameter_list::pair_parameter_type;

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    // [[forcefield.global]]
    // interaciton = "Pair"
    // potential = "LennardJones"
    // table.A.A = {sigma = 1.0, epsilon = 2.0}
    // table.A.B = {sigma = 1.0, epsilon = 2.0} # B.A will be the same
    // table.B.B = {sigma = 1.0, epsilon = 2.0}
    // parameters = [
    //     {index = 0, name = "A"},
    //     {index = 1, name = "A"},
    //     {index = 2, name = "B"},
    //     {index = 3, name = "B"},
    // ]
    std::unordered_map<std::string, pair_parameter_type> table;
    for(const auto& kv : toml::find<toml::table>(global, "table"))
    {
        const auto& p1 = kv.first;
        for(const auto& kv2 : toml::get<toml::table>(kv.second))
        {
            const auto& p2  = kv2.first;
            const auto para = parameter_reader(kv2.second, env);

            const auto key = p1 + std::string(":") + p2;
            table[key] = para;
            if(p1 != p2)
            {
                const auto key_opposite = p2 + std::string(":") + p1;
                if(table.count(key_opposite) != 0)
                {
                    MJOLNIR_LOG_WARN("parameter table does not distinguish "
                                     "two parameters, A.B and B.A.");
                }
                table[key_opposite] = para;
            }
        }
    }

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx  = find_parameter   <std::size_t >(param, env, "index") +
                          find_parameter_or<std::int64_t>(param, env, "offset", 0);
        const auto name = toml::find<std::string>(param, "name");

        params.emplace_back(idx, name);
        MJOLNIR_LOG_INFO("idx = ", idx, ", name = ", name);
    }
    check_parameter_overlap(env, ps, params);

    return make_unique<parameter_list>(std::move(table), std::move(params),
            read_ignore_particles_within(global),
            read_ignored_molecule(global), read_ignored_group(global));
}

template<typename traitsT>
ParameterList<traitsT, ExcludedVolumePotential<typename traitsT::real_type>>
read_excluded_volume_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type      = typename traitsT::real_type;
    using potential_type = ExcludedVolumePotential<real_type>;
    using parameter_list = ExcludedVolumeParameterList<traitsT>;
    using parameter_type = typename parameter_list::parameter_type;

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    const real_type eps = toml::find<real_type>(global, "epsilon");
    MJOLNIR_LOG_INFO("epsilon = ", eps);

    const real_type cutoff = toml::find_or<real_type>(global, "cutoff",
            potential_type::default_cutoff());
    potential_type::set_cutoff_ratio(cutoff);
    MJOLNIR_LOG_INFO("relative cutoff = ", cutoff);

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx    = find_parameter<std::size_t>(param, env, "index") +
                            find_parameter_or<std::int64_t>(param, env, "offset", 0);
        const auto radius = find_parameter<real_type  >(param, env, "radius");

        params.emplace_back(idx, radius);
        MJOLNIR_LOG_INFO("idx = ", idx, ", radius = ", radius);
    }
    check_parameter_overlap(env, ps, params);

    return ParameterList<traitsT, ExcludedVolumePotential<real_type>>(
        make_unique<parameter_list>(eps, cutoff, std::move(params),
            read_ignore_particles_within(global),
            read_ignored_molecule(global), read_ignored_group(global)));
}

template<typename traitsT>
ParameterList<traitsT, InversePowerPotential<typename traitsT::real_type>>
read_inverse_power_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type      = typename traitsT::real_type;
    using potential_type = InversePowerPotential<real_type>;
    using parameter_list = InversePowerParameterList<traitsT>;
    using parameter_type = typename parameter_list::parameter_type;
    using integer_type   = typename potential_type::integer_type;

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    const real_type eps = toml::find<real_type>(global, "epsilon");
    MJOLNIR_LOG_INFO("epsilon = ", eps);

    const integer_type n = toml::find<integer_type>(global, "n");
    MJOLNIR_LOG_INFO("n = ", n);

    const real_type cutoff = toml::find_or<real_type>(global, "cutoff",
            potential_type::default_cutoff(n));
    potential_type::set_cutoff_ratio(cutoff, n);
    MJOLNIR_LOG_INFO("relative cutoff = ", cutoff);

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx    = find_parameter<std::size_t>(param, env, "index") +
                            find_parameter_or<std::int64_t>(param, env, "offset", 0);
        const auto radius = find_parameter<real_type  >(param, env, "radius");

        params.emplace_back(idx, radius);
        MJOLNIR_LOG_INFO("idx = ", idx, ", radius = ", radius);
    }
    check_parameter_overlap(env, ps, params);

    return ParameterList<traitsT, potential_type>(make_unique<parameter_list>(
            eps, n, cutoff, std::move(params),
            read_ignore_particles_within(global),
            read_ignored_molecule(global), read_ignored_group(global)));
}

template<typename traitsT>
ParameterList<traitsT, HardCoreExcludedVolumePotential<typename traitsT::real_type>>
read_hard_core_excluded_volume_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    using real_type      = typename traitsT::real_type;
    using potential_type = HardCoreExcludedVolumePotential<real_type>;
    using parameter_list = HardCoreExcludedVolumeParameterList<traitsT>;
    using parameter_type = typename parameter_list::parameter_type;

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    const real_type eps = toml::find<real_type>(global, "epsilon");
    MJOLNIR_LOG_INFO("epsilon = ", eps);

    const real_type cutoff = toml::find_or<real_type>(global, "cutoff",
            potential_type::default_cutoff());
    potential_type::set_cutoff_ratio(cutoff);
    MJOLNIR_LOG_INFO("relative cutoff = ", cutoff);

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx = find_parameter<std::size_t>(param, env, "index") +
                         find_parameter_or<std::int64_t>(param, env, "offset", 0);

        const auto core_radius =
            find_parameter<real_type>(param, env, "core_radius");
        const auto soft_shell_thickness =
            find_parameter<real_type>(param, env, "soft_shell_thickness");

        params.emplace_back(idx, parameter_type{soft_shell_thickness, core_radius});
        MJOLNIR_LOG_INFO("idx = ", idx, ", core_radius = ", core_radius,
                         ", soft_shell_thickness = ", soft_shell_thickness);
    }
    check_parameter_overlap(env, ps, params);

    return ParameterList<traitsT, potential_type>(make_unique<parameter_list>(
            eps, cutoff, std::move(params), read_ignore_particles_within(global),
            read_ignored_molecule(global), read_ignored_group(global)));
}

template<typename traitsT>
ParameterList<traitsT, LennardJonesPotential<typename traitsT::real_type>>
read_lennard_jones_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    using real_type      = typename traitsT::real_type;
    using potential_type = LennardJonesPotential<real_type>;

    const real_type cutoff = toml::find_or<real_type>(global, "cutoff",
            potential_type::default_cutoff());
    potential_type::set_cutoff_ratio(cutoff);
    MJOLNIR_LOG_INFO("relative cutoff = ", cutoff);

    // check if `table` exists
    if(global.contains("table"))
    {
        return ParameterList<traitsT, potential_type>(
            read_parameter_table<traitsT, potential_type>(global,
                [](const toml::value& val, const toml::value& env) {
                    return std::make_pair(
                        find_parameter<real_type>(val, env, "sigma"),
                        find_parameter<real_type>(val, env, "epsilon"));
                }));
    }
    else
    {
        using parameter_list = LorentzBerthelotRule<traitsT, potential_type>;
        using parameter_type = typename parameter_list::parameter_type;
        // [[forcefield.global]]
        // interaciton = "Pair"
        // potential = "LennardJones"
        // parameters = [
        //     {index = 0, sigma = 1.0, epsilon = 0.2},
        //     # ...
        // ]
        const auto& ps = toml::find<toml::array>(global, "parameters");
        MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

        const auto& env = global.contains("env") ? global.at("env") : toml::value{};

        std::vector<std::pair<std::size_t, parameter_type>> params;
        params.reserve(ps.size());
        for(const auto& param : ps)
        {
            const auto idx     = find_parameter<std::size_t>(param, env, "index") +
                                 find_parameter_or<std::int64_t>(param, env, "offset", 0);
            const auto sigma   = find_parameter<real_type>(param, env, "sigma",   u8"σ");
            const auto epsilon = find_parameter<real_type>(param, env, "epsilon", u8"ε");

            params.emplace_back(idx, parameter_type{sigma, epsilon});
            MJOLNIR_LOG_INFO("idx = ", idx, ", sigma = ", sigma, ", epsilon = ", epsilon);
        }
        check_parameter_overlap(env, ps, params);

        return ParameterList<traitsT, potential_type>(make_unique<parameter_list>(
                std::move(params), read_ignore_particles_within(global),
                read_ignored_molecule(global), read_ignored_group(global)));
    }
}


template<typename traitsT>
ParameterList<traitsT, WCAPotential<typename traitsT::real_type>>
read_wca_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    using real_type      = typename traitsT::real_type;
    using potential_type = WCAPotential<real_type>;

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    if(global.contains("cutoff"))
    {
        MJOLNIR_LOG_WARN("WCA has a exact cutoff distance. While the simulation,"
            " the exact cutoff is used and specified cutoff will be ignored.");
    }

    // check if `table` exists
    if(global.contains("table"))
    {
        return ParameterList<traitsT, potential_type>(
            read_parameter_table<traitsT, potential_type>(global,
                [](const toml::value& val, const toml::value& env) {
                    return std::make_pair(
                        find_parameter<real_type>(val, env, "sigma"),
                        find_parameter<real_type>(val, env, "epsilon"));
                }));
    }
    else
    {
        using parameter_list = LorentzBerthelotRule<traitsT, potential_type>;
        using parameter_type = typename parameter_list::parameter_type;

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

        return ParameterList<traitsT, potential_type>(
            make_unique<parameter_list>(std::move(params),
                read_ignore_particles_within(global),
                read_ignored_molecule(global), read_ignored_group(global)));
    }
}

template<typename traitsT>
ParameterList<traitsT, LennardJonesAttractivePotential<typename traitsT::real_type>>
read_lennard_jones_attractive_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    using real_type      = typename traitsT::real_type;
    using potential_type = LennardJonesAttractivePotential<real_type>;

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    const real_type cutoff = toml::find_or<real_type>(global, "cutoff",
            potential_type::default_cutoff());
    potential_type::set_cutoff_ratio(cutoff);
    MJOLNIR_LOG_INFO("relative cutoff = ", cutoff);

    // check if `table` exists
    if(global.contains("table"))
    {
        return ParameterList<traitsT, potential_type>(
            read_parameter_table<traitsT, potential_type>(global,
                [](const toml::value& val, const toml::value& env) {
                    return std::make_pair(
                        find_parameter<real_type>(val, env, "sigma"),
                        find_parameter<real_type>(val, env, "epsilon"));
                }));
    }
    else
    {
        using parameter_list = LorentzBerthelotRule<traitsT, potential_type>;
        using parameter_type = typename parameter_list::parameter_type;

        const auto& ps = toml::find<toml::array>(global, "parameters");
        MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

        std::vector<std::pair<std::size_t, parameter_type>> params;
        params.reserve(ps.size());
        for(const auto& param : ps)
        {
            const auto idx     = find_parameter<std::size_t>(param, env, "index") +
                                 find_parameter_or<std::int64_t>(param, env, "offset", 0);
            const auto sigma   = find_parameter<real_type>(param, env, "sigma",   u8"σ");
            const auto epsilon = find_parameter<real_type>(param, env, "epsilon", u8"ε");

            params.emplace_back(idx, parameter_type{sigma, epsilon});
            MJOLNIR_LOG_INFO("idx = ", idx, ", sigma = ", sigma, ", epsilon = ", epsilon);
        }

        check_parameter_overlap(env, ps, params);

        return ParameterList<traitsT, potential_type>(make_unique<parameter_list>(
                std::move(params), read_ignore_particles_within(global),
                read_ignored_molecule(global), read_ignored_group(global)));
    }
}

template<typename traitsT>
ParameterList<traitsT, DebyeHuckelPotential<typename traitsT::real_type>>
read_debye_huckel_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    using real_type      = typename traitsT::real_type;
    using potential_type = DebyeHuckelPotential<real_type>;
    using parameter_list = DebyeHuckelParameterList<traitsT>;
    using parameter_type = typename parameter_list::parameter_type;

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    const real_type cutoff = toml::find_or<real_type>(global, "cutoff",
            potential_type::default_cutoff());
//     potential_type::set_cutoff_ratio(cutoff); // will be done in ParameterList
    MJOLNIR_LOG_INFO("relative cutoff = ", cutoff);

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx = find_parameter<std::size_t>(param, env, "index") +
                         find_parameter_or<std::int64_t>(param, env, "offset", 0);
        const auto charge = find_parameter<real_type>(param, env, "charge");

        params.emplace_back(idx, parameter_type{charge});
        MJOLNIR_LOG_INFO("idx = ", idx, ", charge = ", charge);
    }

    check_parameter_overlap(env, ps, params);

    return ParameterList<traitsT, DebyeHuckelPotential<real_type>>(make_unique<parameter_list>(
            cutoff, std::move(params), read_ignore_particles_within(global),
            read_ignored_molecule(global), read_ignored_group(global)));
}

template<typename traitsT>
ParameterList<traitsT, ThreeSPN2ExcludedVolumePotential<typename traitsT::real_type>>
read_3spn2_excluded_volume_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    using real_type      = typename traitsT::real_type;
    using potential_type = ThreeSPN2ExcludedVolumePotential<real_type>;
    using bead_kind      = parameter_3SPN2::bead_kind;
    using parameter_list = ParameterList<traitsT, potential_type>;
    using parameter_type = // here we need to use the specific type of parameter
        typename ThreeSPN2ExcludedVolumeParameterList<traitsT>::parameter_type;

    // cutoff is fixed. it becomes exactly to zero at the cutoff length.

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx = find_parameter<std::size_t>(param, env, "index") +
                         find_parameter_or<std::int64_t>(param, env, "offset", 0);

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
    check_parameter_overlap(env, ps, params);

    return parameter_list(make_unique<ThreeSPN2ExcludedVolumeParameterList<traitsT>>(
            ThreeSPN2ExcludedVolumePotentialParameter<real_type>{},
            params, read_ignore_particles_within(global),
            read_ignored_molecule(global), read_ignored_group(global)));
}

template<typename traitsT>
ParameterList<traitsT, iSoLFAttractivePotential<typename traitsT::real_type>>
read_isolf_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    using real_type      = typename traitsT::real_type;
    using potential_type = iSoLFAttractivePotential<real_type>;
    using parameter_list = iSoLFAttractiveParameterList<traitsT>;
    using parameter_type = typename parameter_list::parameter_type;

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    if(global.contains("cutoff"))
    {
        MJOLNIR_LOG_WARN("iSoLF has an exact cutoff distance. While the simulation,"
            " the exact cutoff is used and specified cutoff will be ignored.");
    }

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        // offset can be a negative value
        const auto idx     = find_parameter<std::int64_t>(param, env, "index") +
                             find_parameter_or<std::int64_t>(param, env, "offset", 0);
        const auto sigma   = find_parameter<real_type>(param, env, "sigma"  );
        const auto epsilon = find_parameter<real_type>(param, env, "epsilon");
        const auto omega   = find_parameter<real_type>(param, env, "omega"  );

        params.emplace_back(idx, parameter_type{sigma, epsilon, omega});
        MJOLNIR_LOG_INFO("idx = ", idx, ", sigma = ", sigma,
                       ", epsilon = ", epsilon, ", omega = ", omega);
    }
    check_parameter_overlap(env, ps, params);

    return ParameterList<traitsT, potential_type>(make_unique<parameter_list>(
            std::move(params), read_ignore_particles_within(global),
            read_ignored_molecule(global), read_ignored_group(global)));

}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , ExcludedVolumePotential<double>> read_excluded_volume_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , ExcludedVolumePotential<float >> read_excluded_volume_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>> read_excluded_volume_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float >> read_excluded_volume_potential(const toml::value&);

extern template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , InversePowerPotential<double>> read_inverse_power_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , InversePowerPotential<float >> read_inverse_power_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, InversePowerPotential<double>> read_inverse_power_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, InversePowerPotential<float >> read_inverse_power_potential(const toml::value&);

extern template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , HardCoreExcludedVolumePotential<double>> read_hard_core_excluded_volume_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , HardCoreExcludedVolumePotential<float >> read_hard_core_excluded_volume_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, HardCoreExcludedVolumePotential<double>> read_hard_core_excluded_volume_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, HardCoreExcludedVolumePotential<float >> read_hard_core_excluded_volume_potential(const toml::value&);

extern template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , LennardJonesPotential<double>> read_lennard_jones_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , LennardJonesPotential<float >> read_lennard_jones_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>> read_lennard_jones_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float >> read_lennard_jones_potential(const toml::value&);

extern template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , LennardJonesAttractivePotential<double>> read_lennard_jones_attractive_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , LennardJonesAttractivePotential<float >> read_lennard_jones_attractive_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesAttractivePotential<double>> read_lennard_jones_attractive_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesAttractivePotential<float >> read_lennard_jones_attractive_potential(const toml::value&);

extern template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , WCAPotential<double>> read_wca_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , WCAPotential<float >> read_wca_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, WCAPotential<double>> read_wca_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, WCAPotential<float >> read_wca_potential(const toml::value&);

extern template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , DebyeHuckelPotential<double>> read_debye_huckel_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , DebyeHuckelPotential<float >> read_debye_huckel_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>> read_debye_huckel_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float >> read_debye_huckel_potential(const toml::value&);

extern template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , ThreeSPN2ExcludedVolumePotential<double>> read_3spn2_excluded_volume_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , ThreeSPN2ExcludedVolumePotential<float >> read_3spn2_excluded_volume_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, ThreeSPN2ExcludedVolumePotential<double>> read_3spn2_excluded_volume_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ThreeSPN2ExcludedVolumePotential<float >> read_3spn2_excluded_volume_potential(const toml::value&);

extern template ParameterList<SimulatorTraits<double, UnlimitedBoundary>       , iSoLFAttractivePotential<double>> read_isolf_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  UnlimitedBoundary>       , iSoLFAttractivePotential<float >> read_isolf_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>, iSoLFAttractivePotential<double>> read_isolf_potential(const toml::value&);
extern template ParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, iSoLFAttractivePotential<float >> read_isolf_potential(const toml::value&);
#endif

} // mjolnir
#endif // MJOLNIR_READ_GLOBAL_POTENTIAL_HPP
