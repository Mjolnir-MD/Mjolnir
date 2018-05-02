#ifndef MJOLNIR_READ_SPATIAL_PARTITION
#define MJOLNIR_READ_SPATIAL_PARTITION
#include <mjolnir/core/UnlimitedGridCellList.hpp>
#include <mjolnir/core/PeriodicGridCellList.hpp>
#include <mjolnir/core/NaivePairCalculation.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/get_toml_value.hpp>

namespace mjolnir
{

template<typename boundaryT, typename traitsT>
struct celllist_dispatcher;

template<typename realT, typename coordT, typename traitsT>
struct celllist_dispatcher<UnlimitedBoundary<realT, coordT>, traitsT>
{
    typedef UnlimitedGridCellList<traitsT> type;
    typedef realT real_type;

    static UnlimitedGridCellList<traitsT> invoke(const real_type mergin)
    {
        return UnlimitedGridCellList<traitsT>(mergin);
    }
};

template<typename realT, typename coordT, typename traitsT>
struct celllist_dispatcher<CubicPeriodicBoundary<realT, coordT>, traitsT>
{
    typedef PeriodicGridCellList<traitsT> type;
    typedef realT real_type;

    static PeriodicGridCellList<traitsT> invoke(const real_type mergin)
    {
        return PeriodicGridCellList<traitsT>(mergin);
    }
};

template<typename traitsT, typename potentialT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_spatial_partition_for_distance(const toml::Table& global, potentialT pot)
{
    typedef typename traitsT::real_type real_type;

    const auto& sp = toml_value_at(
            global, "spatial_partition", "[forcefield.global]"
            ).cast<toml::value_t::Table>();
    const auto  type = toml::get<std::string>(
            toml_value_at(sp, "type", "[forcefield.global]"));

    if(type == "CellList")
    {
        using boundary_type = typename traitsT::boundary_type;
        using dispatcher    = celllist_dispatcher<boundary_type, traitsT>;
        using celllist_type = typename dispatcher::type;

        const auto mg = toml::get<real_type>(
                toml_value_at(sp, "mergin", "[forcefield.global]"));

        return make_unique<GlobalDistanceInteraction<
            traitsT, potentialT, celllist_type>>(
                std::move(pot), dispatcher::invoke(mg));
    }
    else if(type == "VerletList")
    {
        const auto mergin = toml::get<real_type>(toml_value_at(
                    sp, "mergin", "[forcefield.global]"));
        return make_unique<GlobalDistanceInteraction<
            traitsT, potentialT, VerletList<traitsT>>>(
                std::move(pot), VerletList<traitsT>(mergin));
    }
    else if(type == "Naive")
    {
        return make_unique<GlobalDistanceInteraction<
            traitsT, potentialT, NaivePairCalculation<traitsT>>
                >(std::move(pot), NaivePairCalculation<traitsT>());
    }
    else
    {
        throw std::runtime_error("invalid spatial partition type: " + type);
    }
}

} // mjolnir
#endif// MJOLNIR_READ_SPATIAL_PARTITION
