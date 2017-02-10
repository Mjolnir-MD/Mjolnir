#ifndef MJOLNIR_IO_TOML_READ_GLOBAL_INTERACTION
#define MJOLNIR_IO_TOML_READ_GLOBAL_INTERACTION
#include <mjolnir/core/GlobalDistanceInteraction.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <mjolnir/core/UnlimitedGridCellList.hpp>
#include <mjolnir/core/PeriodicGridCellList.hpp>
#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_distance_interaction(const std::string& boundary,
        const typename traitsT::real_type cutoff, const toml::Table& potent)
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

    if(boundary == "Unlimited")
    {
//     auto space = make_unique<VerletList<traitsT>>(cutoff, cutoff * mergin);
        auto space = make_unique<UnlimitedGridCellList<traitsT>>(
                cutoff, cutoff * mergin);
        space->set_except(excepts);
        return make_unique<GlobalDistanceInteraction<traitsT>>(std::move(space));
    }
    else if(boundary == "Periodic")
    {
//         auto space = make_unique<VerletList<traitsT,
//              PeriodicBoundaryXYZ<traitsT>>>(cutoff, cutoff * mergin);
        auto space = make_unique<PeriodicGridCellList<traitsT>>(
                cutoff, cutoff * mergin);
        space->set_except(excepts);
        std::cerr << "periodic exv" << std::endl;
        return make_unique<GlobalDistanceInteraction<traitsT,
               PeriodicBoundaryXYZ<traitsT>>>(std::move(space));
    }
    else
        throw std::runtime_error("unknown boundary condition: " + boundary);
}

template<typename traitsT>
std::unique_ptr<GlobalInteractionBase<traitsT>>
read_global_interaction(const std::string& name,
        const std::unique_ptr<GlobalPotentialBase<traitsT>>& pot,
        const toml::Table& potent)
{
    MJOLNIR_SET_LOGGER("read_toml_file");
    MJOLNIR_LOG_DEBUG("read_global_interaction CALLED");
    MJOLNIR_LOG_INFO("global interaction name", name);

    const std::string boundary = toml::get<toml::String>(potent.at("boundary"));

    if(name == "Global")
    {
        const typename traitsT::real_type cutoff = pot->max_cutoff_length();
        MJOLNIR_LOG_INFO("potential cutoff length", cutoff);

        return read_global_distance_interaction<traitsT>(boundary, cutoff, potent);
    }
    else
        throw std::runtime_error("unknown interaction: " + name);
}


} // mjolnir
#endif /* MJOLNIR_IO_TOML_READ_GLOBAL_INTERACTION */
