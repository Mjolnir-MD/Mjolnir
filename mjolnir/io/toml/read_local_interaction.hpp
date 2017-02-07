#ifndef MJOLNIR_IO_TOML_READ_LOCAL_INTERACTION
#define MJOLNIR_IO_TOML_READ_LOCAL_INTERACTION
#include <mjolnir/core/BondLengthInteraction.hpp>
#include <mjolnir/core/BondAngleInteraction.hpp>
#include <mjolnir/core/DihedralAngleInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>
#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT, 2>>
read_2body_interaction(const std::string& name, const std::string& bound)
{
    if(name == "BondLength")
    {
        if(bound == "Unlimited")
            return make_unique<BondLengthInteraction<traitsT>>();
        else if(bound == "Periodic")
            return make_unique<BondLengthInteraction<traitsT,
                                   PeriodicBoundaryXYZ<traitsT>>>();
        else
            throw std::runtime_error("unknown boundary: " + bound);
    }
    else
    {
        throw std::runtime_error(std::string("unknown interaction: ") + name);
    }
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT, 3>>
read_3body_interaction(const std::string& name, const std::string& bound)
{
    if(name == "BondAngle")
    {
        if(bound == "Unlimited")
            return make_unique<BondAngleInteraction<traitsT>>();
        else if(bound == "Periodic")
            return make_unique<BondAngleInteraction<traitsT,
                                   PeriodicBoundaryXYZ<traitsT>>>();
        else
            throw std::runtime_error("unknown boundary: " + bound);
    }
    else
    {
        throw std::runtime_error(std::string("unknown interaction: ") + name);
    }
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT, 4>>
read_4body_interaction(const std::string& name, const std::string& bound)
{
    if(name == "DihedralAngle")
    {
        if(bound == "Unlimited")
            return make_unique<DihedralAngleInteraction<traitsT>>();
        else if(bound == "Periodic")
            return make_unique<DihedralAngleInteraction<traitsT,
                                   PeriodicBoundaryXYZ<traitsT>>>();
        else
            throw std::runtime_error("unknown boundary: " + bound);
    }
    else
    {
        throw std::runtime_error(std::string("unknown interaction: ") + name);
    }
}

} // mjolnir
#endif /* MJOLNIR_IO_TOML_READ_LOCAL_INTERACTION */
