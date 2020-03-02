#include <mjolnir/interaction/local/DihedralAngleInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// ClementiDihedral
template class DihedralAngleInteraction<SimulatorTraits<double, UnlimitedBoundary>, ClementiDihedralPotential<double>>;
template class DihedralAngleInteraction<SimulatorTraits<float,  UnlimitedBoundary>, ClementiDihedralPotential<float> >;
template class DihedralAngleInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ClementiDihedralPotential<double>>;
template class DihedralAngleInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ClementiDihedralPotential<float> >;

// gaussian
template class DihedralAngleInteraction<SimulatorTraits<double, UnlimitedBoundary>, PeriodicGaussianPotential<double>>;
template class DihedralAngleInteraction<SimulatorTraits<float,  UnlimitedBoundary>, PeriodicGaussianPotential<float> >;
template class DihedralAngleInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, PeriodicGaussianPotential<double>>;
template class DihedralAngleInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, PeriodicGaussianPotential<float> >;

// FLP dihedral
template class DihedralAngleInteraction<SimulatorTraits<double, UnlimitedBoundary>, FlexibleLocalDihedralPotential<double>>;
template class DihedralAngleInteraction<SimulatorTraits<float,  UnlimitedBoundary>, FlexibleLocalDihedralPotential<float> >;
template class DihedralAngleInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, FlexibleLocalDihedralPotential<double>>;
template class DihedralAngleInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, FlexibleLocalDihedralPotential<float> >;

// cosine
template class DihedralAngleInteraction<SimulatorTraits<double, UnlimitedBoundary>, CosinePotential<double>>;
template class DihedralAngleInteraction<SimulatorTraits<float,  UnlimitedBoundary>, CosinePotential<float> >;
template class DihedralAngleInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>>;
template class DihedralAngleInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> >;

} // mjolnir
