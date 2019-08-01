#ifndef MJOLNIR_INPUT_READ_SPATIAL_PARTITION_HPP
#define MJOLNIR_INPUT_READ_SPATIAL_PARTITION_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/GlobalInteractionBase.hpp>
#include <mjolnir/interaction/global/GlobalPairInteraction.hpp>

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
//
// XXX note for specialist:
// After constructing neighbor-list, we can calculate forcefield parameter, e.g.
// multiply of charges of each pair, in order to accelerate force calculation.
// But to store those parameters in an efficient way, we need to pass the type
// of parameters to the cell-list. Here, parameterT represents the type.
template<typename boundaryT, typename traitsT, typename parameterT>
struct celllist_dispatcher;

// implementation for Unlimited Boundary case.
template<typename realT, typename coordT, typename traitsT, typename parameterT>
struct celllist_dispatcher<
    UnlimitedBoundary<realT, coordT>, traitsT, parameterT
    >
{
    using real_type = realT;
    using type      = UnlimitedGridCellList<traitsT, parameterT>;

    static UnlimitedGridCellList<traitsT, parameterT>
    invoke(const real_type margin)
    {
        return UnlimitedGridCellList<traitsT, parameterT>(margin);
    }
};

// implementation for Cuboidal Periodic Boundary case.
template<typename realT, typename coordT, typename traitsT, typename parameterT>
struct celllist_dispatcher<
    CuboidalPeriodicBoundary<realT, coordT>, traitsT, parameterT
    >
{
    using real_type = realT;
    using type      = PeriodicGridCellList<traitsT, parameterT>;

    static PeriodicGridCellList<traitsT, parameterT>
    invoke(const real_type margin)
    {
        return PeriodicGridCellList<traitsT, parameterT>(margin);
    }
};

// ---------------------------------------------------------------------------
// It reads spatial partition that is dedicated for a GlobalPotential.
// In most of the cases, "type" would be a "CellList".
//
// TODO: this assumes the global interaction is GlobalPairInteraction.
//       After some other interactions are added, we would need to change
//       the interface.
template<typename traitsT, typename potentialT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_spatial_partition(const toml::value& global, potentialT&& pot)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type      = typename traitsT::real_type;
    using parameter_type = typename potentialT::pair_parameter_type;

    const auto& sp   = toml::find<toml::value>(global, "spatial_partition");
    const auto  type = toml::find<std::string>(sp,     "type");

    if(type == "CellList")
    {
        using boundary_type = typename traitsT::boundary_type;
        using dispatcher    = celllist_dispatcher<boundary_type, traitsT, parameter_type>;
        using celllist_type = typename dispatcher::type;

        const auto margin = toml::find<real_type>(sp, "margin");
        MJOLNIR_LOG_NOTICE("-- Spatial Partition is CellList "
                           "with relative margin = ", margin);
        return make_unique<
            GlobalPairInteraction<traitsT, potentialT, celllist_type>
            >(std::forward<potentialT>(pot), dispatcher::invoke(margin));
    }
    else if(type == "VerletList")
    {
        using verlet_list_type = VerletList<traitsT, parameter_type>;

        const auto margin = toml::find<real_type>(sp, "margin");
        MJOLNIR_LOG_NOTICE("-- Spatial Partition is VerletList "
                           "with relative margin = ", margin);
        return make_unique<
            GlobalPairInteraction<traitsT, potentialT, verlet_list_type>
            >(std::forward<potentialT>(pot), verlet_list_type(margin));
    }
    else if(type == "Naive")
    {
        MJOLNIR_LOG_NOTICE("-- No Spatial Partition. Calculate all the possible pairs.");
        using naive_pair_type = NaivePairCalculation<traitsT, parameter_type>;

        return make_unique<
            GlobalPairInteraction<traitsT, potentialT, naive_pair_type>
            >(std::forward<potentialT>(pot), naive_pair_type());
    }
    else
    {
        throw std::runtime_error(toml::format_error("[error] "
            "mjolnir::read_spatial_partition: unknown option appeared",
            toml::find<toml::value>(sp, "type"),
            "expected \"CellList\" or \"VerletList\""));
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
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, DebyeHuckelPotential<double>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, DebyeHuckelPotential<double>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, DebyeHuckelPotential<double>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, DebyeHuckelPotential<double>&&);

extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, DebyeHuckelPotential<float>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, DebyeHuckelPotential<float>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, DebyeHuckelPotential<float>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, DebyeHuckelPotential<float>&&);

// EXV
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, ExcludedVolumePotential<double>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, ExcludedVolumePotential<double>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, ExcludedVolumePotential<double>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, ExcludedVolumePotential<double>&&);

extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, ExcludedVolumePotential<float>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, ExcludedVolumePotential<float>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, ExcludedVolumePotential<float>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, ExcludedVolumePotential<float>&&);

// L-J
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, LennardJonesPotential<double>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, LennardJonesPotential<double>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, LennardJonesPotential<double>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, LennardJonesPotential<double>&&);

extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, LennardJonesPotential<float>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, LennardJonesPotential<float>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, LennardJonesPotential<float>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, LennardJonesPotential<float>&&);

// UL-J
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, UniformLennardJonesPotential<double>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, UniformLennardJonesPotential<double>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, UniformLennardJonesPotential<double>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, UniformLennardJonesPotential<double>&&);

extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, UniformLennardJonesPotential<float>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_spatial_partition(const toml::value&, UniformLennardJonesPotential<float>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, UniformLennardJonesPotential<float>&&);
extern template std::unique_ptr<GlobalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_spatial_partition(const toml::value&, UniformLennardJonesPotential<float>&&);

}
#endif

#endif// MJOLNIR_READ_SPATIAL_PARTITION
