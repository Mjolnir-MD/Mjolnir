#ifndef MJOLNIR_IO_TOML_READ_SPATIAL_PARTITION
#define MJOLNIR_IO_TOML_READ_SPATIAL_PARTITION
#include <mjolnir/core/UnlimitedGridCellList.hpp>
#include <mjolnir/core/PeriodicGridCellList.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT, typename spaceT, typename boundaryT>
struct read_spatial_partition_impl;

template<typename traitsT, typename spaceT, typename boundaryT>
spaceT
read_spatial_partition(const toml::Table& pot,
        const typename traitsT::real_type cutoff)
{
    return read_spatial_partition_impl<traitsT, spaceT, boundaryT>::invoke(
            pot, cutoff);
}

template<typename traitsT, typename boundaryT>
struct read_spatial_partition_impl<traitsT,
    VerletList<traitsT, boundaryT>, boundaryT>
{
    typedef typename traitsT::real_type real_type;

    static VerletList<traitsT, boundaryT>
    invoke(const toml::Table& potent, const real_type cutoff)
    {
        MJOLNIR_SET_LOGGER("read_toml_file");
        typename traitsT::real_type mergin = 1.0;
        try{mergin = toml::get<toml::Float>(potent.at("mergin"));}
        catch(std::exception& except){mergin = 1.0;}
        MJOLNIR_LOG_INFO("mergin rate", mergin);

        std::vector<std::vector<std::size_t>> excepts;
        auto elist = toml::get<toml::Array<toml::Array<toml::Integer>>>(
                potent.at("excepts"));
        MJOLNIR_LOG_INFO("except list size", elist.size());

        for(auto iter = elist.cbegin(); iter != elist.cend(); ++iter)
        {
            std::vector<std::size_t> l; l.reserve(iter->size());
            for(auto j = iter->cbegin(); j != iter->cend(); ++j)
                l.push_back(static_cast<std::size_t>(*j));

            excepts.emplace_back(l);
        }
        auto retval = VerletList<traitsT, boundaryT>(cutoff, cutoff * mergin);
        retval.set_except(excepts);
        return retval;
    }
};

template<typename traitsT, typename boundaryT>
struct read_spatial_partition_impl<traitsT,
    UnlimitedGridCellList<traitsT, boundaryT>, boundaryT>
{
    typedef typename traitsT::real_type real_type;

    static UnlimitedGridCellList<traitsT, boundaryT>
    invoke(const toml::Table& potent, const real_type cutoff)
    {
        MJOLNIR_SET_LOGGER("read_toml_file");
        typename traitsT::real_type mergin = 1.0;
        try{mergin = toml::get<toml::Float>(potent.at("mergin"));}
        catch(std::exception& except){mergin = 1.0;}
        MJOLNIR_LOG_INFO("mergin rate", mergin);

        std::vector<std::vector<std::size_t>> excepts;
        auto elist = toml::get<toml::Array<toml::Array<toml::Integer>>>(
                potent.at("excepts"));
        MJOLNIR_LOG_INFO("except list size", elist.size());

        for(auto iter = elist.cbegin(); iter != elist.cend(); ++iter)
        {
            std::vector<std::size_t> l; l.reserve(iter->size());
            for(auto j = iter->cbegin(); j != iter->cend(); ++j)
                l.push_back(static_cast<std::size_t>(*j));

            excepts.emplace_back(l);
        }
        auto retval = UnlimitedGridCellList<traitsT>(cutoff, cutoff * mergin);
        retval.set_except(excepts);
        return retval;
    }
};

template<typename traitsT, typename boundaryT>
struct read_spatial_partition_impl<traitsT,
    PeriodicGridCellList<traitsT, boundaryT>, boundaryT>
{
    typedef typename traitsT::real_type real_type;

    static PeriodicGridCellList<traitsT, boundaryT>
    invoke(const toml::Table& potent, const real_type cutoff)
    {
        MJOLNIR_SET_LOGGER("read_toml_file");
        typename traitsT::real_type mergin = 1.0;
        try{mergin = toml::get<toml::Float>(potent.at("mergin"));}
        catch(std::exception& except){mergin = 1.0;}
        MJOLNIR_LOG_INFO("mergin rate", mergin);

        std::vector<std::vector<std::size_t>> excepts;
        auto elist = toml::get<toml::Array<toml::Array<toml::Integer>>>(
                potent.at("excepts"));
        MJOLNIR_LOG_INFO("except list size", elist.size());

        for(auto iter = elist.cbegin(); iter != elist.cend(); ++iter)
        {
            std::vector<std::size_t> l; l.reserve(iter->size());
            for(auto j = iter->cbegin(); j != iter->cend(); ++j)
                l.push_back(static_cast<std::size_t>(*j));

            excepts.emplace_back(l);
        }
        auto retval = PeriodicGridCellList<traitsT>(cutoff, cutoff * mergin);
        retval.set_except(excepts);
        return retval;
    }
};

}//mjolnir
#endif/* MJOLNIR_IO_TOML_READ_SPATIAL_PARTITION*/
