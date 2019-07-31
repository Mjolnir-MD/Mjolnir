#ifndef MJOLNIR_INPUT_READ_GLOBAL_POTENTIAL_HPP
#define MJOLNIR_INPUT_READ_GLOBAL_POTENTIAL_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/input/utility.hpp>
#include <mjolnir/potential/global/ExcludedVolumePotential.hpp>
#include <mjolnir/potential/global/LennardJonesPotential.hpp>
#include <mjolnir/potential/global/UniformLennardJonesPotential.hpp>
#include <mjolnir/potential/global/DebyeHuckelPotential.hpp>
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

inline IgnoreMolecule<typename Topology::molecule_id_type>
read_ignored_molecule(const toml::value& ignore)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

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

inline IgnoreGroup<typename Topology::group_id_type>
read_ignored_group(const toml::value& ignore)
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

    if(ignore.as_table().count("group") == 0)
    {
        MJOLNIR_LOG_NOTICE("No `ignore.group` is provided. "
                           "All the groups are taken into account.");
        assert(ignores.empty());
        return IgnoreGroup<group_id_type>(ignores);
    }

    const auto group = toml::find(ignore, "group");
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
        const auto inter = toml::find(group, "inter");
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

template<typename realT>
ExcludedVolumePotential<realT>
read_excluded_volume_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = realT;

    const auto& ignore = toml::find<toml::value>(global, "ignore");

    const auto ignore_particle_within = toml::find<
        std::map<std::string, std::size_t>>(ignore, "particles_within");
    for(const auto& connection : ignore_particle_within)
    {
        MJOLNIR_LOG_INFO("particles that have connection ", connection.first,
            " within ", connection.second, " will be ignored");
    }

    const auto& env = global.as_table().count("env") == 1 ?
                      global.as_table().at("env") : toml::value{};

    const real_type eps = toml::find<real_type>(global, "epsilon");
    MJOLNIR_LOG_INFO("epsilon = ", eps);

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    using parameter_type = typename ExcludedVolumePotential<realT>::parameter_type;

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx    = find_parameter<std::size_t>(param, env, "index");
        const auto radius = find_parameter<real_type  >(param, env, "radius");

        params.emplace_back(idx, radius);
        MJOLNIR_LOG_INFO("idx = ", idx, ", radius = ", radius);
    }

    return ExcludedVolumePotential<realT>(
        eps, params, ignore_particle_within,
        read_ignored_molecule(ignore), read_ignored_group(ignore));
}

template<typename realT>
LennardJonesPotential<realT>
read_lennard_jones_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = realT;

    const auto& ignore = toml::find<toml::value>(global, "ignore");

    const auto ignore_particle_within = toml::find<
        std::map<std::string, std::size_t>>(ignore, "particles_within");
    for(const auto& connection : ignore_particle_within)
    {
        MJOLNIR_LOG_INFO("particles that have connection ", connection.first,
            " within ", connection.second, " will be ignored");
    }

    const auto& env = global.as_table().count("env") == 1 ?
                      global.as_table().at("env") : toml::value{};

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    using parameter_type = typename LennardJonesPotential<realT>::parameter_type;

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx     = find_parameter<std::size_t>(param, env, "index");
        const auto sigma   = find_parameter<real_type>(param, env, "sigma",   u8"σ");
        const auto epsilon = find_parameter<real_type>(param, env, "epsilon", u8"ε");

        params.emplace_back(idx, parameter_type{sigma, epsilon});
        MJOLNIR_LOG_INFO("idx = ", idx, ", sigma = ", sigma, ", epsilon = ", epsilon);
    }

    return LennardJonesPotential<realT>(
        std::move(params), ignore_particle_within,
        read_ignored_molecule(ignore), read_ignored_group(ignore));
}

template<typename realT>
UniformLennardJonesPotential<realT>
read_uniform_lennard_jones_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = realT;

    const auto& ignore = toml::find<toml::value>(global, "ignore");

    const auto ignore_particle_within = toml::find<
        std::map<std::string, std::size_t>>(ignore, "particles_within");
    for(const auto& connection : ignore_particle_within)
    {
        MJOLNIR_LOG_INFO("particles that have connection ", connection.first,
            " within ", connection.second, " will be ignored");
    }

    const auto& env = global.as_table().count("env") == 1 ?
                      global.as_table().at("env") : toml::value{};

    const auto sigma   = toml::expect<real_type>(global, u8"σ").or_other(
                         toml::expect<real_type>(global, "sigma")).unwrap();
    const auto epsilon = toml::expect<real_type>(global, u8"ε").or_other(
                         toml::expect<real_type>(global, "epsilon")).unwrap();

    MJOLNIR_LOG_INFO("sigma   = ", sigma);
    MJOLNIR_LOG_INFO("epsilon = ", epsilon);

    using parameter_type = typename UniformLennardJonesPotential<realT>::parameter_type;
    std::vector<std::pair<std::size_t, parameter_type>> params;
    if(global.as_table().count("parameters") == 1)
    {
        for(const auto& param : toml::find<toml::array>(global, "parameters"))
        {
            const auto idx = find_parameter<std::size_t>(param, env, "index");
            params.emplace_back(idx, parameter_type{});
        }
    }
    return UniformLennardJonesPotential<realT>(
        sigma, epsilon, params, ignore_particle_within,
        read_ignored_molecule(ignore), read_ignored_group(ignore));
}

template<typename realT>
DebyeHuckelPotential<realT>
read_debye_huckel_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = realT;

    const auto& ignore = toml::find<toml::value>(global, "ignore");

    const auto ignore_particle_within = toml::find<
        std::map<std::string, std::size_t>>(ignore, "particles_within");
    for(const auto& connection : ignore_particle_within)
    {
        MJOLNIR_LOG_INFO("particles that have connection ", connection.first,
            " within ", connection.second, " will be ignored");
    }

    const auto& env = global.as_table().count("env") == 1 ?
                      global.as_table().at("env") : toml::value{};

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    using parameter_type = typename DebyeHuckelPotential<realT>::parameter_type;

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx    = find_parameter<std::size_t>(param, env, "index");
        const auto charge = find_parameter<real_type  >(param, env, "charge");

        params.emplace_back(idx, parameter_type{charge});
        MJOLNIR_LOG_INFO("idx = ", idx, ", charge = ", charge);
    }
    return DebyeHuckelPotential<realT>(
        std::move(params), ignore_particle_within,
        read_ignored_molecule(ignore), read_ignored_group(ignore));
}

template<typename realT>
ThreeSPN2ExcludedVolumePotential<realT>
read_3spn2_excluded_volume_potential(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using parameter_type =
        typename ThreeSPN2ExcludedVolumePotential<realT>::parameter_type;
    using bead_kind = parameter_3SPN2::bead_kind;

    const auto& env = global.as_table().count("env") == 1 ?
                      global.as_table().at("env") : toml::value{};

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::pair<std::size_t, parameter_type>> params;
    params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto idx  = find_parameter<std::size_t>(param, env, "index");
        const auto kind = find_parameter<std::string>(param, env, "kind");

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

    IgnoreGroup<typename Topology::group_id_type> ignore_grp({});
    if(global.as_table().count("ignore") == 1)
    {
        const auto& ignore = toml::find<toml::value>(global, "ignore");
        ignore_grp = read_ignored_group(ignore);
    }
    return ThreeSPN2ExcludedVolumePotential<realT>(params, ignore_grp);
}


#ifdef MJOLNIR_SEPARATE_BUILD
extern template ExcludedVolumePotential<double> read_excluded_volume_potential(const toml::value& global);
extern template ExcludedVolumePotential<float > read_excluded_volume_potential(const toml::value& global);

extern template LennardJonesPotential<double> read_lennard_jones_potential(const toml::value& global);
extern template LennardJonesPotential<float > read_lennard_jones_potential(const toml::value& global);

extern template UniformLennardJonesPotential<double> read_uniform_lennard_jones_potential(const toml::value& global);
extern template UniformLennardJonesPotential<float > read_uniform_lennard_jones_potential(const toml::value& global);

extern template DebyeHuckelPotential<double> read_debye_huckel_potential(const toml::value& global);
extern template DebyeHuckelPotential<float > read_debye_huckel_potential(const toml::value& global);
#endif

} // mjolnir
#endif // MJOLNIR_READ_GLOBAL_POTENTIAL_HPP
