#include <mjolnir/forcefield/iSoLF/iSoLFAttractivePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class iSoLFAttractiveParameterList<SimulatorTraits<double, UnlimitedBoundary>       >;
template class iSoLFAttractiveParameterList<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class iSoLFAttractiveParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class iSoLFAttractiveParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
