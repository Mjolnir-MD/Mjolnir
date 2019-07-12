#include <mjolnir/core/NaivePairCalculation.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>, empty_t>;
template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>, empty_t>;
template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, empty_t>;
template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, empty_t>;

template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>, double>;
template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>, float >;
template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, double>;
template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, float >;

template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>, std::pair<double, double>>;
template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>, std::pair<float , float >>;
template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, std::pair<double, double>>;
template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, std::pair<float , float >>;
} // mjolnir
