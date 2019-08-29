#ifndef MJOLNIR_INPUT_READ_SPATIAL_PARTITION_HPP
#define MJOLNIR_INPUT_READ_SPATIAL_PARTITION_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/SpatialPartitionBase.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/UnlimitedGridCellList.hpp>
#include <mjolnir/core/PeriodicGridCellList.hpp>
#include <mjolnir/core/NaivePairCalculation.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// ---------------------------------------------------------------------------
// It constructs cell list depending on the boundary condition.
// Since the most efficient implementation changes depending on the boundary
// condition, there is `celllist_dispatcher` that creates different cell list
// depending on the boundary condition.

template<typename boundaryT>
struct celllist_dispatcher;

// implementation for Unlimited Boundary case.
template<typename realT, typename coordT>
struct celllist_dispatcher<UnlimitedBoundary<realT, coordT>>
{
    using real_type = realT;

    template<typename traitsT, typename potentialT>
    static std::unique_ptr<UnlimitedGridCellList<traitsT, potentialT>>
    invoke(const real_type margin)
    {
        return make_unique<UnlimitedGridCellList<traitsT, potentialT>>(margin);
    }
};

// implementation for Cuboidal Periodic Boundary case.
template<typename realT, typename coordT>
struct celllist_dispatcher<CuboidalPeriodicBoundary<realT, coordT>>
{
    using real_type = realT;

    template<typename traitsT, typename potentialT>
    static std::unique_ptr<PeriodicGridCellList<traitsT, potentialT>>
    invoke(const real_type margin)
    {
        return make_unique<PeriodicGridCellList<traitsT, potentialT>>(margin);
    }
};

// ---------------------------------------------------------------------------
// It reads spatial partition that is dedicated for a GlobalPotential.
// In most of the cases, "type" would be a "CellList".

template<typename traitsT, typename potentialT>
SpatialPartition<traitsT, potentialT>
read_spatial_partition(const toml::value& global)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = typename traitsT::real_type;

    const auto& sp   = toml::find<toml::value>(global, "spatial_partition");
    const auto  type = toml::find<std::string>(sp,     "type");

    if(type == "CellList")
    {
        using boundary_type = typename traitsT::boundary_type;

        const auto margin = toml::find<real_type>(sp, "margin");
        MJOLNIR_LOG_NOTICE("-- Spatial Partition is CellList "
                           "with relative margin = ", margin);

        return SpatialPartition<traitsT, potentialT>(
                celllist_dispatcher<boundary_type>::template
                invoke<traitsT, potentialT>(margin));
    }
    else if(type == "VerletList")
    {
        using verlet_list_type = VerletList<traitsT, potentialT>;

        const auto margin = toml::find<real_type>(sp, "margin");
        MJOLNIR_LOG_NOTICE("-- Spatial Partition is VerletList "
                           "with relative margin = ", margin);

        return SpatialPartition<traitsT, potentialT>(
                make_unique<verlet_list_type>(margin));
    }
    else if(type == "Naive")
    {
        MJOLNIR_LOG_NOTICE("-- No Spatial Partition. "
                           "Calculate all the possible pairs.");

        return SpatialPartition<traitsT, potentialT>(
                make_unique<NaivePairCalculation<traitsT, potentialT>>());
    }
    else
    {
        throw std::runtime_error(toml::format_error("[error] "
            "mjolnir::read_spatial_partition: unknown option appeared",
            toml::find(sp, "type"), "expected \"CellList\" or \"VerletList\""));
    }
}

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD

#include <mjolnir/core/SimulatorTraits.hpp>

#include <mjolnir/potential/global/DebyeHuckelPotential.hpp>
#include <mjolnir/potential/global/ExcludedVolumePotential.hpp>
#include <mjolnir/potential/global/LennardJonesPotential.hpp>
#include <mjolnir/potential/global/UniformLennardJonesPotential.hpp>

namespace mjolnir
{
// D-H
extern template SpatialPartition<SimulatorTraits<double, UnlimitedBoundary>       , DebyeHuckelPotential<double>        > read_spatial_partition(const toml::value&);
extern template SpatialPartition<SimulatorTraits<float,  UnlimitedBoundary>       , DebyeHuckelPotential<float >        > read_spatial_partition(const toml::value&);
extern template SpatialPartition<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>        > read_spatial_partition(const toml::value&);
extern template SpatialPartition<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float >        > read_spatial_partition(const toml::value&);

// EXV
extern template SpatialPartition<SimulatorTraits<double, UnlimitedBoundary>       , ExcludedVolumePotential<double>     > read_spatial_partition(const toml::value&);
extern template SpatialPartition<SimulatorTraits<float,  UnlimitedBoundary>       , ExcludedVolumePotential<float >     > read_spatial_partition(const toml::value&);
extern template SpatialPartition<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>     > read_spatial_partition(const toml::value&);
extern template SpatialPartition<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float >     > read_spatial_partition(const toml::value&);

// L-J
extern template SpatialPartition<SimulatorTraits<double, UnlimitedBoundary>       , LennardJonesPotential<double>       > read_spatial_partition(const toml::value&);
extern template SpatialPartition<SimulatorTraits<float,  UnlimitedBoundary>       , LennardJonesPotential<float >       > read_spatial_partition(const toml::value&);
extern template SpatialPartition<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>       > read_spatial_partition(const toml::value&);
extern template SpatialPartition<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float >       > read_spatial_partition(const toml::value&);

// UL-J
extern template SpatialPartition<SimulatorTraits<double, UnlimitedBoundary>       , UniformLennardJonesPotential<double>> read_spatial_partition(const toml::value&);
extern template SpatialPartition<SimulatorTraits<float,  UnlimitedBoundary>       , UniformLennardJonesPotential<float >> read_spatial_partition(const toml::value&);
extern template SpatialPartition<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>> read_spatial_partition(const toml::value&);
extern template SpatialPartition<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float >> read_spatial_partition(const toml::value&);

}
#endif

#endif// MJOLNIR_READ_SPATIAL_PARTITION
