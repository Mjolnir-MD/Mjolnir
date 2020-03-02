#include <mjolnir/interaction/local/DirectionalContactInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, CosinePotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , CosinePotential<float> , GaussianPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, CosinePotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , CosinePotential<float> , GaussianPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, GaussianPotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , GaussianPotential<float> , GaussianPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, GaussianPotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , GaussianPotential<float> , GaussianPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, UniformPotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , UniformPotential<float> , GaussianPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, UniformPotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , UniformPotential<float> , GaussianPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, GaussianPotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , GaussianPotential<float> , GaussianPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, GaussianPotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , GaussianPotential<float> , GaussianPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, UniformPotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , UniformPotential<float> , GaussianPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, UniformPotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , UniformPotential<float> , GaussianPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, CosinePotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , CosinePotential<float> , GaussianPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, CosinePotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , CosinePotential<float> , GaussianPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, UniformPotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , UniformPotential<float> , GaussianPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, UniformPotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , UniformPotential<float> , GaussianPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, CosinePotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , CosinePotential<float> , GaussianPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, CosinePotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , CosinePotential<float> , GaussianPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, GaussianPotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , GaussianPotential<float> , GaussianPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, GaussianPotential<double>, GaussianPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , GaussianPotential<float> , GaussianPotential<float> >;

// ---------------------------------------------------------

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, CosinePotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , CosinePotential<float> , GoContactPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, CosinePotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , CosinePotential<float> , GoContactPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, GaussianPotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , GaussianPotential<float> , GoContactPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, GaussianPotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , GaussianPotential<float> , GoContactPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, UniformPotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , UniformPotential<float> , GoContactPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, UniformPotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , UniformPotential<float> , GoContactPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, GaussianPotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , GaussianPotential<float> , GoContactPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, GaussianPotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , GaussianPotential<float> , GoContactPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, UniformPotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , UniformPotential<float> , GoContactPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, UniformPotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , UniformPotential<float> , GoContactPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, CosinePotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , CosinePotential<float> , GoContactPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, CosinePotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , CosinePotential<float> , GoContactPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, UniformPotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , UniformPotential<float> , GoContactPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, UniformPotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , UniformPotential<float> , GoContactPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, CosinePotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , CosinePotential<float> , GoContactPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, CosinePotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , CosinePotential<float> , GoContactPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, GaussianPotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , GaussianPotential<float> , GoContactPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, GaussianPotential<double>, GoContactPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , GaussianPotential<float> , GoContactPotential<float> >;

// --------------------------------------------

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, CosinePotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , CosinePotential<float> , UniformPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, CosinePotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , CosinePotential<float> , UniformPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, GaussianPotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , GaussianPotential<float> , UniformPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, GaussianPotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , GaussianPotential<float> , UniformPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, UniformPotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , UniformPotential<float> , UniformPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, UniformPotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , UniformPotential<float> , UniformPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, GaussianPotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , GaussianPotential<float> , UniformPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, GaussianPotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , GaussianPotential<float> , UniformPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, UniformPotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , UniformPotential<float> , UniformPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, UniformPotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , UniformPotential<float> , UniformPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, CosinePotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , CosinePotential<float> , UniformPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, CosinePotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , CosinePotential<float> , UniformPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, UniformPotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , UniformPotential<float> , UniformPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, UniformPotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , UniformPotential<float> , UniformPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, CosinePotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , CosinePotential<float> , UniformPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, CosinePotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , CosinePotential<float> , UniformPotential<float> >;

template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, GaussianPotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , GaussianPotential<float> , UniformPotential<float> >;
template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, GaussianPotential<double>, UniformPotential<double>>;
template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , GaussianPotential<float> , UniformPotential<float> >;

} // mjolnir
