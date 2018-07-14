#ifndef MJOLNIR_READ_SPATIAL_PARTITION
#define MJOLNIR_READ_SPATIAL_PARTITION
#include <mjolnir/core/UnlimitedGridCellList.hpp>
#include <mjolnir/core/PeriodicGridCellList.hpp>
#include <mjolnir/core/NaivePairCalculation.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

template<typename boundaryT, typename traitsT>
struct celllist_dispatcher;

template<typename realT, typename coordT, typename traitsT>
struct celllist_dispatcher<UnlimitedBoundary<realT, coordT>, traitsT>
{
    typedef UnlimitedGridCellList<traitsT> type;
    typedef realT real_type;

    static UnlimitedGridCellList<traitsT> invoke(const real_type margin)
    {
        return UnlimitedGridCellList<traitsT>(margin);
    }
};

template<typename realT, typename coordT, typename traitsT>
struct celllist_dispatcher<CuboidalPeriodicBoundary<realT, coordT>, traitsT>
{
    typedef PeriodicGridCellList<traitsT> type;
    typedef realT real_type;

    static PeriodicGridCellList<traitsT> invoke(const real_type margin)
    {
        return PeriodicGridCellList<traitsT>(margin);
    }
};

template<typename traitsT, typename potentialT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_spatial_partition(const toml::Table& global, potentialT pot)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_spatial_partition(), 0);
    typedef typename traitsT::real_type real_type;

    const auto& sp   = get_toml_value<toml::Table>(
            global, "spatial_partition", "[forcefield.global]");
    const auto  type = get_toml_value<std::string>(
            sp, "type", "[forcefield.global]");

    if(type == "CellList")
    {
        MJOLNIR_SCOPE(type == "CellList", 1);
        using boundary_type = typename traitsT::boundary_type;
        using dispatcher    = celllist_dispatcher<boundary_type, traitsT>;
        using celllist_type = typename dispatcher::type;

        const auto mg =
            get_toml_value<real_type>(sp, "margin", "[forcefield.global]");
        MJOLNIR_LOG_INFO("margin = ", mg);

        return make_unique<GlobalPairInteraction<
            traitsT, potentialT, celllist_type>>(
                std::move(pot), dispatcher::invoke(mg));
    }
    else if(type == "VerletList")
    {
        MJOLNIR_SCOPE(type == "VerletList", 1);

        const auto margin = get_toml_value<real_type>(
                    sp, "margin", "[forcefield.global]");
        MJOLNIR_LOG_INFO("margin = ", margin);

        return make_unique<GlobalPairInteraction<
            traitsT, potentialT, VerletList<traitsT>>>(
                std::move(pot), VerletList<traitsT>(margin));
    }
    else if(type == "Naive")
    {
        MJOLNIR_SCOPE(type == "Naive", 1);
        return make_unique<GlobalPairInteraction<
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
