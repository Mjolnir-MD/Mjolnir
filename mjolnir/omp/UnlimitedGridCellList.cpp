#include <mjolnir/omp/UnlimitedGridCellList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, DebyeHuckelPotential<OpenMPSimulatorTraits<double, UnlimitedBoundary>>>;
template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, DebyeHuckelPotential<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>>;

template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, ExcludedVolumePotential<OpenMPSimulatorTraits<double, UnlimitedBoundary>>>;
template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, ExcludedVolumePotential<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>>;

template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, LennardJonesPotential<OpenMPSimulatorTraits<double, UnlimitedBoundary>>>;
template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, LennardJonesPotential<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>>;
} // mjolnir
