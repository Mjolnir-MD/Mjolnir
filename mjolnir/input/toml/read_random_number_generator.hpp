#ifndef MJOLNIR_TOML_READ_RANDOM_NUMBER_GENERATOR
#define MJOLNIR_TOML_READ_RANDOM_NUMBER_GENERATOR
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <mjolnir/util/logger.hpp>
#include <toml/toml.hpp>
#include <memory>

namespace mjolnir
{

template<typename traitsT>
std::shared_ptr<RandomNumberGenerator<traitsT>>
read_random_number_generator(const toml::Table& sim)
{
    MJOLNIR_SET_LOGGER("read_toml_file");
    MJOLNIR_LOG_DEBUG("read_random_number_generator CALLED");
    const std::uint32_t seed = toml::get<toml::Integer>(sim.at("seed"));
    MJOLNIR_LOG_INFO("seed", seed);
    MJOLNIR_LOG_DEBUG("read_random_number_generator RETURNED");
    return std::make_shared<RandomNumberGenerator<traitsT>>(seed);
}

}// mjolnir
#endif /* MJOLNIR_TOML_READ_RANDOM_NUMBER_GENERATOR */
