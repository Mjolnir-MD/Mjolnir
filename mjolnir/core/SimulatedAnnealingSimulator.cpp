#include <mjolnir/core/SimulatedAnnealingSimulator.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
// BAOAB
template class SimulatedAnnealingSimulator<SimulatorTraits<double, UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >, LinearScheduler>;
template class SimulatedAnnealingSimulator<SimulatorTraits<float,  UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >, LinearScheduler>;
template class SimulatedAnnealingSimulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>, LinearScheduler>;
template class SimulatedAnnealingSimulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>, LinearScheduler>;
// Langevin
template class SimulatedAnnealingSimulator<SimulatorTraits<double, UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >, LinearScheduler>;
template class SimulatedAnnealingSimulator<SimulatorTraits<float , UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >, LinearScheduler>;
template class SimulatedAnnealingSimulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>, LinearScheduler>;
template class SimulatedAnnealingSimulator<SimulatorTraits<float , CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>, LinearScheduler>;
} // mjolnir
