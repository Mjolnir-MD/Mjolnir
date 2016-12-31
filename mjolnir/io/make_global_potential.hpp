#ifndef MJOLNIR_MAKE_GLOBAL_POTENTIAL
#define MJOLNIR_MAKE_GLOBAL_POTENTIAL
#include <mjolnir/core/GlobalPotentialBase.hpp>
#include <mjolnir/core/LennardJonesPotential.hpp>
#include <mjolnir/core/ExcludedVolumePotential.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <stdexcept>
#include <string>

namespace mjolnir
{

template<typename traitsT, typename ... Ts>
std::unique_ptr<GlobalPotentialBase<traitsT>>
make_global_potential(const std::string& name, Ts&&... args)
{
    if(name == "ExcludedVolume")
    {
        return mjolnir::make_unique<ExcludedVolumePotential<traitsT>>(
                std::forward<Ts>(args)...);
    }
    else if(name == "LennardJones")
    {
        return mjolnir::make_unique<LennardJonesPotential<traitsT>>(
                std::forward<Ts>(args)...);
    }
    else
    {
        throw std::invalid_argument("Unknown Potential: " + name);
    }
}

}// mjolnir
#endif /* MJOLNIR_MAKE_GLOBAL_POTENTIAL */
