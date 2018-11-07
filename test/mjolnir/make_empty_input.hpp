#ifndef MJOLNIR_TEST_MAKE_EMPTY_INPUT_HPP
#define MJOLNIR_TEST_MAKE_EMPTY_INPUT_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/util/string.hpp>

namespace mjolnir
{
namespace test
{

// make a minimal input file that contains no particle, no forcefield.
// It contains only parameters that are needed to pass read_input() functions.
inline toml::Table make_empty_input()
{
    using namespace mjolnir::literals::string_literals;

    toml::Table general;
    general["precision"]     = "double"_s;
    general["boundary"]      = "Unlimited"_s;
    general["output_prefix"] = "nothing"_s;
    general["output_path"]   = "./"_s;

    toml::Table units;
    units["length"] = "angstrom"_s;
    units["energy"] = "kcal/mol"_s;

    // Mjolnir does not consider running simulation without simulator
    // (that does not make sense!) . Here, temporary set MD simulator with
    // Newtonian dynamics.
    toml::Table simulator;
    simulator["type"]       = "Molecular Dynamics"_s;
    simulator["scheme"]     = "Newtonian"_s;
    simulator["total_step"] = 1;
    simulator["save_step"]  = 1;
    simulator["delta_t"]    = 0.1;

    // Mjolnir requires a system to simulate, even if it has no particle.
    toml::Table system;
    system["attributes"] = toml::Table{{"temperature"_s, toml::value(300.0)}};
    system["boundary"]   = toml::Table(); // no boundary
    system["particles"]  = toml::Array(); // no particle
    toml::Array systems(1, std::move(system));

    // "empty forcefields" completely make sense because essentially any kind of
    // forcefields are a sum of several potential terms. Forcefield that has
    // zero term means an ideal system having no interaction.
    toml::Array forcefields;

    toml::Table input;
    input["general"]     = std::move(general);
    input["units"]       = std::move(units);
    input["simulator"]   = std::move(simulator);
    input["systems"]     = std::move(systems);
    input["forcefields"] = std::move(forcefields);

    return input;
}

} // test
} // mjolnir
#endif// MJOLNIR_TEST_MAKE_EMPTY_INPUT_HPP
