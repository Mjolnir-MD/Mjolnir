#include <mjolnir/omp/UnlimitedGridCellList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, DebyeHuckelPotential<OpenMPSimulatorTraits<double, UnlimitedBoundary>>>;
template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, DebyeHuckelPotential<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>>;

template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, ExcludedVolumePotential<double>>;
template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, ExcludedVolumePotential<float >>;

template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, LennardJonesPotential<double>>;
template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, LennardJonesPotential<float >>;

template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, UniformLennardJonesPotential<double>>;
template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, UniformLennardJonesPotential<float >>;
} // mjolnir
