#ifndef MJOLNIR_READ_SPATIAL_PARTITION
#define MJOLNIR_READ_SPATIAL_PARTITION
#include <mjolnir/core/UnlimitedGridCellList.hpp>
#include <mjolnir/core/PeriodicGridCellList.hpp>
#include <mjolnir/core/NaivePairCalculation.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <mjolnir/util/make_unique.hpp>

namespace mjolnir
{

template<typename boundaryT, typename traitsT>
struct celllist_dispatcher;

template<typename realT, typename coordT, typename traitsT>
struct celllist_dispatcher<UnlimitedBoundary<realT, coordT>, traitsT>
{
    typedef UnlimitedGridCellList<traitsT> type;
    typedef realT real_type;
    static UnlimitedGridCellList<traitsT>
    invoke(const real_type cutoff, const real_type mergin)
    {
        return UnlimitedGridCellList<traitsT>{cutoff, mergin};
    }
};

template<typename realT, typename coordT, typename traitsT>
struct celllist_dispatcher<CubicPeriodicBoundary<realT, coordT>, traitsT>
{
    typedef PeriodicGridCellList<traitsT> type;
    typedef realT real_type;
    static PeriodicGridCellList<traitsT>
    invoke(const real_type cutoff, const real_type mergin)
    {
        return PeriodicGridCellList<traitsT>{cutoff, mergin};
    }
};

template<typename partitionT>
partitionT
read_exception_information(const toml::Table& global, partitionT&& sp)
{
    const auto& params = global.at("parameters").cast<toml::value_t::Array>();
    for(const auto& tab : params)
    {
        const auto& info = tab.cast<toml::value_t::Table>();
        const auto  idx = toml::get<std::size_t>(info.at("index"));
        sp.chain_index(idx) = toml::get<std::size_t>(info.at("chain"));
        for(auto exc : toml::get<std::vector<std::size_t>>(
                    info.at("except_chains")))
        {
            sp.except_chains(idx).push_back(exc);
        }
        for(auto exb : toml::get<std::vector<std::size_t>>(
                    info.at("except_beads")))
        {
            sp.except_indices(idx).push_back(exb);
        }
    }
    return sp;
}


template<typename traitsT, typename potentialT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_spatial_partition_for_distance(const toml::Table& global, potentialT&& pot)
{
    typedef typename traitsT::real_type real_type;

    const auto& sp = global.at("spatial_partition").cast<toml::value_t::Table>();
    const auto  type = toml::get<std::string>(sp.at("type"));
    if(type == "CellList")
    {
        typedef typename traitsT::boundary_type boundary_type;
        typedef typename celllist_dispatcher<boundary_type, traitsT>::type
                celllist_type;

        const auto co = pot.max_cutoff_length();
        const auto mg = toml::get<real_type>(sp.at("mergin"));
        return make_unique<GlobalDistanceInteraction<
            traitsT, potentialT, celllist_type>>(std::move(pot),
                celllist_dispatcher<boundary_type, traitsT>::invoke(co, mg));
    }
    else if(type == "VerletList")
    {
        const auto cutoff = pot.max_cutoff_length();
        const auto mergin = toml::get<real_type>(sp.at("mergin"));
        return make_unique<GlobalDistanceInteraction<
            traitsT, potentialT, VerletList<traitsT>>
            >(std::move(pot), VerletList<traitsT>{cutoff, mergin});
    }
    else if(type == "Naive")
    {
        return make_unique<GlobalDistanceInteraction<
            traitsT, potentialT, NaivePairCalculation<traitsT>>
            >(std::move(pot), NaivePairCalculation<traitsT>{});
    }
    else
    {
        throw std::runtime_error("invalid spatial partition type: " + type);
    }
}

template<typename traitsT, typename potentialT>
std::unique_ptr<GlobalInteractionBase<traitT>>
read_spatial_partition_for_externalMU(const toml::Table& global, potentialT&& pot)
{
    return make_unique<GlobalExternalInteraction<
      traitsT, potentialT, SpatialPartitionForMU<traitsT>>
		       >(std::move(pot), SpatialPartitionForMU<traitsT>{});
}
}
#endif// MJOLNIR_READ_SPATIAL_PARTITION
