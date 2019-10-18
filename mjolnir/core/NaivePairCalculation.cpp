#include <mjolnir/core/NaivePairCalculation.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>       , DebyeHuckelPotential<SimulatorTraits<double, UnlimitedBoundary>       >>;
template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>       , DebyeHuckelPotential<SimulatorTraits<float,  UnlimitedBoundary>       >>;
template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>       , ExcludedVolumePotential<double>>;
template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>       , ExcludedVolumePotential<float >>;
template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>>;
template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float >>;

template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>       , LennardJonesPotential<double>>;
template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>       , LennardJonesPotential<float >>;
template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>>;
template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float >>;

template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>       , UniformLennardJonesPotential<double>>;
template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>       , UniformLennardJonesPotential<float >>;
template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>>;
template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float >>;
} // mjolnir
