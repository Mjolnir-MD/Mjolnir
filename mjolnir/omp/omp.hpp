#ifndef MJOLNIR_OMP_OMP_HPP
#define MJOLNIR_OMP_OMP_HPP

// This file is a meta-header file that just includes everything that are needed
// to use OpenMP implementation. This file is introduced in order to make it
// simple to turn on/off OpenMP support by preprocessor macro, like in the
// following way.
//
// ```cpp
// #ifdef WITH_OPENMP
// #include <mjolnir/omp/omp.hpp>
// #endif
// ```
//
// By doing this, the general, standalone implementation does not need to
// consider what files should be included to use OpenMP.

#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/omp/RandomNumberGenerator.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/BondLengthInteraction.hpp>
#include <mjolnir/omp/BondLengthGoContactInteraction.hpp>
#include <mjolnir/omp/ContactInteraction.hpp>
#include <mjolnir/omp/GoContactInteraction.hpp>
#include <mjolnir/omp/BondAngleInteraction.hpp>
#include <mjolnir/omp/DihedralAngleInteraction.hpp>
#include <mjolnir/omp/ThreeSPN2BaseStackingInteraction.hpp>
#include <mjolnir/omp/ProteinDNANonSpecificInteraction.hpp>
#include <mjolnir/omp/GlobalPairInteraction.hpp>
#include <mjolnir/omp/GlobalPairExcludedVolumeInteraction.hpp>
#include <mjolnir/omp/GlobalPairLennardJonesInteraction.hpp>
#include <mjolnir/omp/ThreeSPN2BasePairInteraction.hpp>
#include <mjolnir/omp/ThreeSPN2CrossStackingInteraction.hpp>
#include <mjolnir/omp/PWMcosInteraction.hpp>
#include <mjolnir/omp/PositionRestraintInteraction.hpp>
#include <mjolnir/omp/ExternalDistanceInteraction.hpp>
#include <mjolnir/omp/UnlimitedGridCellList.hpp>
#include <mjolnir/omp/PeriodicGridCellList.hpp>
#include <mjolnir/omp/UnderdampedLangevinIntegrator.hpp>
#include <mjolnir/omp/GFWNPTLangevinIntegrator.hpp>
#include <mjolnir/omp/BAOABLangevinIntegrator.hpp>
#include <mjolnir/omp/gBAOABLangevinIntegrator.hpp>
#include <mjolnir/omp/GJFNVTLangevinIntegrator.hpp>
#include <mjolnir/omp/SystemMotionRemover.hpp>

#endif// MJOLNIR_OMP_OMP_HPP
