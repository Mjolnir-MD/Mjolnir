#include <mjolnir/core/VerletList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class VerletList<SimulatorTraits<double, UnlimitedBoundary>,        DebyeHuckelPotential<double>, true>;
template class VerletList<SimulatorTraits<float,  UnlimitedBoundary>,        DebyeHuckelPotential<float >, true>;
template class VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>, true>;
template class VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float >, true>;

template class VerletList<SimulatorTraits<double, UnlimitedBoundary>,        ExcludedVolumePotential<double>, true>;
template class VerletList<SimulatorTraits<float,  UnlimitedBoundary>,        ExcludedVolumePotential<float >, true>;
template class VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>, true>;
template class VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float >, true>;

template class VerletList<SimulatorTraits<double, UnlimitedBoundary>,        LennardJonesPotential<double>, true>;
template class VerletList<SimulatorTraits<float,  UnlimitedBoundary>,        LennardJonesPotential<float >, true>;
template class VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>, true>;
template class VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float >, true>;

template class VerletList<SimulatorTraits<double, UnlimitedBoundary>,        UniformLennardJonesPotential<double>, true>;
template class VerletList<SimulatorTraits<float,  UnlimitedBoundary>,        UniformLennardJonesPotential<float >, true>;
template class VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>, true>;
template class VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float >, true>;
} // mjolnir
