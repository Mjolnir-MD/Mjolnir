#include <mjolnir/forcefield/external/RectangularBoxInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class RectangularBoxInteraction<SimulatorTraits<double, UnlimitedBoundary>,        ExcludedVolumeWallPotential<double>>;
template class RectangularBoxInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        ExcludedVolumeWallPotential<float >>;
template class RectangularBoxInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<double>>;
template class RectangularBoxInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<float >>;

template class RectangularBoxInteraction<SimulatorTraits<double, UnlimitedBoundary>,        LennardJonesWallPotential<double>>;
template class RectangularBoxInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        LennardJonesWallPotential<float >>;
template class RectangularBoxInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesWallPotential<double>>;
template class RectangularBoxInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesWallPotential<float >>;

} // mjolnir
