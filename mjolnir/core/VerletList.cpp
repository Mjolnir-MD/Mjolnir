#include <mjolnir/core/VerletList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class VerletList<SimulatorTraits<double, UnlimitedBoundary>,        DebyeHuckelPotential<double>>;
template class VerletList<SimulatorTraits<float,  UnlimitedBoundary>,        DebyeHuckelPotential<float >>;
template class VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>>;
template class VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float >>;

template class VerletList<SimulatorTraits<double, UnlimitedBoundary>,        ExcludedVolumePotential<double>>;
template class VerletList<SimulatorTraits<float,  UnlimitedBoundary>,        ExcludedVolumePotential<float >>;
template class VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>>;
template class VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float >>;

template class VerletList<SimulatorTraits<double, UnlimitedBoundary>,        LennardJonesPotential<double>>;
template class VerletList<SimulatorTraits<float,  UnlimitedBoundary>,        LennardJonesPotential<float >>;
template class VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>>;
template class VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float >>;

template class VerletList<SimulatorTraits<double, UnlimitedBoundary>,        UniformLennardJonesPotential<double>>;
template class VerletList<SimulatorTraits<float,  UnlimitedBoundary>,        UniformLennardJonesPotential<float >>;
template class VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>>;
template class VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float >>;
} // mjolnir
