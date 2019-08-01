#include <mjolnir/core/VerletList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class VerletList<SimulatorTraits<double, UnlimitedBoundary>, empty_t>;
template class VerletList<SimulatorTraits<float,  UnlimitedBoundary>, empty_t>;
template class VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>, empty_t>;
template class VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, empty_t>;

template class VerletList<SimulatorTraits<double, UnlimitedBoundary>, double>;
template class VerletList<SimulatorTraits<float,  UnlimitedBoundary>, float >;
template class VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>, double>;
template class VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, float >;

template class VerletList<SimulatorTraits<double, UnlimitedBoundary>, std::pair<double, double>>;
template class VerletList<SimulatorTraits<float,  UnlimitedBoundary>, std::pair<float , float >>;
template class VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>, std::pair<double, double>>;
template class VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, std::pair<float , float >>;
} // mjolnir
