#ifndef MJOLNIR_INPUT_READ_RANDOM_NUMBER_GENERATOR_HPP
#define MJOLNIR_INPUT_READ_RANDOM_NUMBER_GENERATOR_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <mjolnir/core/MsgPackLoader.hpp>

namespace mjolnir
{

template<typename traitsT>
RandomNumberGenerator<traitsT> read_rng(const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    std::uint32_t seed = 0;
    if(simulator.as_table().count("integrator") != 0)
    {
        const auto& integrator = toml::find(simulator, "integrator");
        if(integrator.as_table().count("seed") != 0)
        {
            MJOLNIR_LOG_WARN("deprecated: put `seed` under [simulator] table.");
            MJOLNIR_LOG_WARN("deprecated: ```toml");
            MJOLNIR_LOG_WARN("deprecated: [simulator]");
            MJOLNIR_LOG_WARN("deprecated: seed = 12345");
            MJOLNIR_LOG_WARN("deprecated: ```");
            seed = toml::find<std::uint32_t>(integrator, "seed");
        }
        else
        {
            seed = toml::find<std::uint32_t>(simulator, "seed");
        }
    }
    else
    {
        if(simulator.at("seed").is_integer())
        {
            seed = toml::find<std::uint32_t>(simulator, "seed");
        }
        else // load from saved checkpoint file
        {
            const std::string fname = toml::find<std::string>(simulator, "seed");
            MsgPackLoader<traitsT> loader;
            MJOLNIR_LOG_NOTICE("RNG is loaded from ", fname);
            return loader.load_rng(fname);
        }
    }
    MJOLNIR_LOG_NOTICE("seed is ", seed);
    return RandomNumberGenerator<traitsT>(seed);
}

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template RandomNumberGenerator<SimulatorTraits<double, UnlimitedBoundary>       > read_rng(const toml::value&);
extern template RandomNumberGenerator<SimulatorTraits<float,  UnlimitedBoundary>       > read_rng(const toml::value&);
extern template RandomNumberGenerator<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_rng(const toml::value&);
extern template RandomNumberGenerator<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_rng(const toml::value&);
}
#endif // MJOLNIR_SEPARATE_BUILD

#endif// MJOLNIR_INPUT_READ_RANDOM_NUMBER_GENERATOR_HPP
