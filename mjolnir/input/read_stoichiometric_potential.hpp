#ifndef MJOLNIR_INPUT_READ_STOICHIOMETRIC_POTENTIAL_HPP
#define MJOLNIR_INPUT_READ_STOICHIOMETRIC_POTENTIAL_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/input/utility.hpp>
#include <mjolnir/input/read_global_potential.hpp>
#include <mjolnir/forcefield/stoichiometric/ParameterList.hpp>
#include <mjolnir/forcefield/stoichiometric/StoichiometricUniformCubicPanPotential.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// ---------------------------------------------------------------------------
// stoichiometric interaction
// ---------------------------------------------------------------------------

template<typename traitsT>
std::pair<StoichiometricUniformCubicPanPotential<typename traitsT::real_type>,
    StoichiometricEmptyCombinationRule<
        traitsT, StoichiometricUniformCubicPanPotential<typename traitsT::real_type>
        >>
read_stoichiometric_uniform_cubic_pan_potential(const toml::value& global,
        const std::string& particle_a_name, const std::string& particle_b_name)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type           = typename traitsT::real_type;
    using potential_type      = StoichiometricUniformCubicPanPotential<real_type>;
    using parameter_list_type = StoichiometricEmptyCombinationRule<traitsT, potential_type>;

    const auto& env = global.contains("env") ? global.at("env") : toml::value{};

    const real_type v0    = toml::find<real_type>(global, "v0");
    MJOLNIR_LOG_INFO("v0 = ", v0);

    const real_type range = toml::find<real_type>(global, "range");
    MJOLNIR_LOG_INFO("range = ", range);

    potential_type potential(v0, range);

    const auto& ps = toml::find<toml::array>(global, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<std::size_t> a_indices;
    std::vector<std::size_t> b_indices;
    for(const auto& param : ps)
    {
        const auto idx  = find_parameter<std::size_t>(param, env, "index") +
                          find_parameter_or<std::int64_t>(param, env, "offset", 0);
        const auto name = toml::find<std::string>(param, "name");

        if(name == particle_a_name)
        {
            a_indices.push_back(idx);
        }
        else if(name == particle_b_name)
        {
            b_indices.push_back(idx);
        }
        else
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_stoichiometric_interaction: unknown particle kind ",
                toml::find<toml::value>(global, "kinds"), "here expected value is \"",
                {particle_a_name, "\" or \"", particle_b_name, "\"."}));
        }
        MJOLNIR_LOG_INFO("idx = ", idx, ", name = ", name);
    }
    return std::make_pair(std::move(potential),
        parameter_list_type(std::move(a_indices), std::move(b_indices),
            read_ignore_particles_within(global),
            read_ignored_molecule(global), read_ignored_group(global)));
}

} // mjolnir

#endif //  MJOLNIR_INPUT_READ_STOICHIOMETRIC_POTENTIAL_HPP

