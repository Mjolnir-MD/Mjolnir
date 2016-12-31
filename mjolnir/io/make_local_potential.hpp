#ifndef MJOLNIR_MAKE_LOCAL_POTENTIAL
#define MJOLNIR_MAKE_LOCAL_POTENTIAL
#include <mjolnir/core/LocalPotentialBase.hpp>
#include <mjolnir/core/HarmonicPotential.hpp>
#include <mjolnir/core/ClementiDihedralPotential.hpp>
#include <mjolnir/core/Go1012ContactPotential.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <stdexcept>
#include <string>

namespace mjolnir
{

template<typename traitsT, typename ... Ts>
std::unique_ptr<LocalPotentialBase<traitsT>>
make_local_potential(const std::string& name, Ts&&... args)
{
    if(name == "Harmonic")
    {
        return mjolnir::make_unique<HarmonicPotential<traitsT>>(
                std::forward<Ts>(args)...);
    }
    else if(name == "ClementiDihedral")
    {
        return mjolnir::make_unique<ClementiDihedralPotential<traitsT>>(
                std::forward<Ts>(args)...);
    }
    else if(name == "Go1012Contact")
    {
        return mjolnir::make_unique<Go1012ContactPotential<traitsT>>(
                std::forward<Ts>(args)...);
    }
    else
    {
        throw std::invalid_argument("Unknown Potential: " + name);
    }
}


}// mjolnir
#endif /* MJOLNIR_MAKE_LOCAL_POTENTIAL */
