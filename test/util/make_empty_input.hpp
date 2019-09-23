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
inline toml::value make_empty_input()
{
    using namespace toml::literals::toml_literals;
    return u8R"(
[files]
output.prefix = "empty"
output.format = "xyz"
output.path   = "./"

[units]
length = "angstrom"
energy = "kcal/mol"

[simulator]
type            = "MolecularDynamics"
precision       = "double"
boundary_type   = "Unlimited"
seed            = 1
total_step      = 1
save_step       = 1
delta_t         = 0.1
integrator.type = "VelocityVerlet"

[[systems]]
attributes.temperature = 300.0
boundary_shape         = {}
particles              = []

[[forcefields]]
# nothing!
)"_toml;
}

} // test
} // mjolnir
#endif// MJOLNIR_TEST_MAKE_EMPTY_INPUT_HPP
