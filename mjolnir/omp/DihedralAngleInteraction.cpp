#include <mjolnir/omp/DihedralAngleInteraction.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{

// ClementiDihedral
template class DihedralAngleInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>, ClementiDihedralPotential<double>>;
template class DihedralAngleInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, ClementiDihedralPotential<float> >;
template class DihedralAngleInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, ClementiDihedralPotential<double>>;
template class DihedralAngleInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, ClementiDihedralPotential<float> >;

// gaussian
template class DihedralAngleInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>, PeriodicGaussianPotential<double>>;
template class DihedralAngleInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, PeriodicGaussianPotential<float> >;
template class DihedralAngleInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, PeriodicGaussianPotential<double>>;
template class DihedralAngleInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, PeriodicGaussianPotential<float> >;

// FLP dihedral
template class DihedralAngleInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>, FlexibleLocalDihedralPotential<double>>;
template class DihedralAngleInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, FlexibleLocalDihedralPotential<float> >;
template class DihedralAngleInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, FlexibleLocalDihedralPotential<double>>;
template class DihedralAngleInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, FlexibleLocalDihedralPotential<float> >;

// cosine
template class DihedralAngleInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>, CosinePotential<double>>;
template class DihedralAngleInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, CosinePotential<float> >;
template class DihedralAngleInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>>;
template class DihedralAngleInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> >;

} // mjolnir
