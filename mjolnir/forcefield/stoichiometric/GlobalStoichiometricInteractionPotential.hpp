#ifndef MJOLNIR_FORCEFIELD_GLOBAL_STOICHIOMETRIC_INTERACTION_POTENTIAL_HPP
#define MJOLNIR_FORCEFIELD_GLOBAL_STOICHIOMETRIC_INTERACTION_POTENTIAL_HPP
#include <mjolnir/core/System.hpp>

namespace mjolnir
{

template<typename traitsT>
class GlobalStoichiometricInteractionPotential
{
  public:
    using traits_type   = traitsT;

    struct pair_parameter_type
    {};

  public:
    GlobalStoichiometricInteractionPotential() = default;
};

} // mjolnir

#endif // MJOLNIR_POTENTIAL_GLOBAL_STOICHIOMETRIC_INTERACTION_POTENTIAL_HPP
