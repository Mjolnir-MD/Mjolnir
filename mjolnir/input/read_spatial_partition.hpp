#ifndef MJOLNIR_READ_SPATIAL_PARTITION
#define MJOLNIR_READ_SPATIAL_PARTITION
#include <mjolnir/core/UnlimitedGridCellList.hpp>
#include <mjolnir/core/PeriodicGridCellList.hpp>
#include <mjolnir/core/NaivePairCalculation.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <mjolnir/core/ImplicitMembraneList.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/input/get_toml_value.hpp>

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
    const auto& params = toml_value_at(global, "parameters", "[forcefield.global]"
            ).cast<toml::value_t::Array>();
    for(const auto& tab : params)
    {
        const auto& info = tab.cast<toml::value_t::Table>();
        const auto  idx = toml::get<std::size_t>(toml_value_at(
                    info, "index", "<anonymous> in parameters"));
        sp.chain_index(idx) = toml::get<std::size_t>(toml_value_at(
                    info, "chain", "<anonymous> in parameters"));

        for(auto exc : toml::get<std::vector<std::size_t>>(
                    toml_value_at(info, "except_chains",
                        "<anonymous> in parameters")))
        {
            sp.except_chains(idx).push_back(exc);
        }
        for(auto exb : toml::get<std::vector<std::size_t>>(
                    toml_value_at(info, "except_beads",
                        "<anonymous> in parameters")))
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

    const auto& sp = toml_value_at(
            global, "spatial_partition", "[forcefield.global]"
            ).cast<toml::value_t::Table>();
    const auto  type = toml::get<std::string>(
            toml_value_at(sp, "type", "[forcefield.global]"));
    if(type == "CellList")
    {
        typedef typename traitsT::boundary_type boundary_type;
        typedef typename celllist_dispatcher<boundary_type, traitsT>::type
                celllist_type;

        const auto co = pot.max_cutoff_length();
        const auto mg = toml::get<real_type>(toml_value_at(
                    sp, "mergin", "[forcefield.global]"));
        return make_unique<GlobalDistanceInteraction<
            traitsT, potentialT, celllist_type>>(std::move(pot),
                read_exception_information(global,
                    celllist_dispatcher<boundary_type, traitsT>::invoke(co, mg)));
    }
    else if(type == "VerletList")
    {
        const auto cutoff = pot.max_cutoff_length();
        const auto mergin = toml::get<real_type>(toml_value_at(
                    sp, "mergin", "[forcefield.global]"));
        return make_unique<GlobalDistanceInteraction<
            traitsT, potentialT, VerletList<traitsT>>>(std::move(pot),
                read_exception_information(
                    global, VerletList<traitsT>{cutoff, mergin}));
    }
    else if(type == "Naive")
    {
        return make_unique<GlobalDistanceInteraction<
            traitsT, potentialT, NaivePairCalculation<traitsT>>
            >(std::move(pot), read_exception_information(
                    global, NaivePairCalculation<traitsT>{}));
    }
    else
    {
        throw std::runtime_error("invalid spatial partition type: " + type);
    }
}

template<typename traitsT, typename potentialT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_spatial_partition_for_implicit_membrane(const toml::Table& global, potentialT&& pot)
{
    const auto cutoff = pot.max_cutoff_length();
    return make_unique<ZaxisExternalForceInteraction<
      traitsT, potentialT, ImplicitMembraneList<traitsT>>
		       >(std::move(pot), ImplicitMembraneList<traitsT>{cutoff});
}
}
#endif// MJOLNIR_READ_SPATIAL_PARTITION
