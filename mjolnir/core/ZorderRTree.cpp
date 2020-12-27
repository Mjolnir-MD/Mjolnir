#include <mjolnir/core/ZorderRTree.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class ZorderRTree<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ZorderRTree<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ZorderRTree<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ZorderRTree<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ZorderRTree<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ZorderRTree<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ZorderRTree<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
template class ZorderRTree<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

template class ZorderRTree<SimulatorTraits<double, UnlimitedBoundary>, DebyeHuckelPotential<SimulatorTraits<double, UnlimitedBoundary>>>;
template class ZorderRTree<SimulatorTraits<float,  UnlimitedBoundary>, DebyeHuckelPotential<SimulatorTraits<float,  UnlimitedBoundary>>>;

template class ZorderRTree<SimulatorTraits<double, UnlimitedBoundary>, ExcludedVolumePotential<SimulatorTraits<double, UnlimitedBoundary>>>;
template class ZorderRTree<SimulatorTraits<float,  UnlimitedBoundary>, ExcludedVolumePotential<SimulatorTraits<float,  UnlimitedBoundary>>>;

template class ZorderRTree<SimulatorTraits<double, UnlimitedBoundary>, LennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>>>;
template class ZorderRTree<SimulatorTraits<float,  UnlimitedBoundary>, LennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>>>;

template class ZorderRTree<SimulatorTraits<double, UnlimitedBoundary>, UniformLennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>>>;
template class ZorderRTree<SimulatorTraits<float,  UnlimitedBoundary>, UniformLennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>>>;
} // mjolnir
