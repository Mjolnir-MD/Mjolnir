#include <mjolnir/forcefield/iSoLF/iSoLFAttractivePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class iSoLFAttractivePotential<SimulatorTraits<double, UnlimitedBoundary>       >;
template class iSoLFAttractivePotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
template class iSoLFAttractivePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
template class iSoLFAttractivePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
