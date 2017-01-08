#ifndef MJOLNIR_TOML_READ_RANDOM_NUMBER_GENERATOR
#define MJOLNIR_TOML_READ_RANDOM_NUMBER_GENERATOR
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <memory>
#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT>
std::shared_ptr<RandomNumberGenerator<traitsT>>
read_random_number_generator(const toml::Table& sim)
{
    const std::uint32_t seed = toml::get<toml::Integer>(sim.at("seed"));
    return std::make_shared<RandomNumberGenerator<traitsT>>(seed);
}

}// mjolnir
#endif /* MJOLNIR_TOML_READ_RANDOM_NUMBER_GENERATOR */
