#include <mjolnir/forcefield/3SPN2/ThreeSPN2BaseBaseInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

template class ThreeSPN2BaseBaseInteraction<SimulatorTraits<double, UnlimitedBoundary>,        UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>,        typename ThreeSPN2BaseBaseInteractionPotential<double>::pair_parameter_type>>;
template class ThreeSPN2BaseBaseInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>,        typename ThreeSPN2BaseBaseInteractionPotential<float>::pair_parameter_type> >;
template class ThreeSPN2BaseBaseInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, PeriodicGridCellList <SimulatorTraits<double, CuboidalPeriodicBoundary>, typename ThreeSPN2BaseBaseInteractionPotential<double>::pair_parameter_type>>;
template class ThreeSPN2BaseBaseInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, PeriodicGridCellList <SimulatorTraits<float,  CuboidalPeriodicBoundary>, typename ThreeSPN2BaseBaseInteractionPotential<float>::pair_parameter_type> >;

template class ThreeSPN2BaseBaseInteraction<SimulatorTraits<double, UnlimitedBoundary>,        VerletList<SimulatorTraits<double, UnlimitedBoundary>,        typename ThreeSPN2BaseBaseInteractionPotential<double>::pair_parameter_type>>;
template class ThreeSPN2BaseBaseInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        VerletList<SimulatorTraits<float,  UnlimitedBoundary>,        typename ThreeSPN2BaseBaseInteractionPotential<float>::pair_parameter_type> >;
template class ThreeSPN2BaseBaseInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>, typename ThreeSPN2BaseBaseInteractionPotential<double>::pair_parameter_type>>;
template class ThreeSPN2BaseBaseInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, typename ThreeSPN2BaseBaseInteractionPotential<float>::pair_parameter_type> >;

template class ThreeSPN2BaseBaseInteraction<SimulatorTraits<double, UnlimitedBoundary>,        NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>,        typename ThreeSPN2BaseBaseInteractionPotential<double>::pair_parameter_type>>;
template class ThreeSPN2BaseBaseInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>,        typename ThreeSPN2BaseBaseInteractionPotential<float>::pair_parameter_type> >;
template class ThreeSPN2BaseBaseInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, typename ThreeSPN2BaseBaseInteractionPotential<double>::pair_parameter_type>>;
template class ThreeSPN2BaseBaseInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, typename ThreeSPN2BaseBaseInteractionPotential<float>::pair_parameter_type> >;

} // mjolnir
