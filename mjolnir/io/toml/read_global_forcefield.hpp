#ifndef MJOLNIR_IO_TOML_READ_GLOBAL_FORCEFIELD
#define MJOLNIR_IO_TOML_READ_GLOBAL_FORCEFIELD
#include <mjolnir/core/GlobalDistanceInteraction.hpp>
#include <mjolnir/core/GlobalForceField.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <mjolnir/core/UnlimitedGridCellList.hpp>
#include <mjolnir/potential/ExcludedVolumePotential.hpp>
#include <mjolnir/potential/LennardJonesPotential.hpp>
#include <mjolnir/util/zip_iterator.hpp>
#include <mjolnir/util/make_zip.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>

#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<GlobalPotentialBase<traitsT>>
read_global_potential(const std::string& name, const toml::Table& potent)
{
    MJOLNIR_SET_LOGGER("read_toml_file");
    MJOLNIR_LOG_DEBUG("read_global_potential CALLED");
    MJOLNIR_LOG_INFO("global potential name", name);

    if(name == "ExcludedVolume")
    {
        toml::Array<toml::Table> param_table =
            toml::get<toml::Array<toml::Table>>(potent.at("parameters"));
        MJOLNIR_LOG_INFO("parameter table size", param_table.size());

        std::vector<typename ExcludedVolumePotential<traitsT>::parameter_type>
            params(param_table.size());
        for(auto iter = make_zip(param_table.cbegin(), params.begin());
                iter != make_zip(param_table.cend(), params.end()); ++iter)
        {
            const typename traitsT::real_type sigma =
                toml::get<toml::Float>(get<0>(iter)->at("sigma"));
            MJOLNIR_LOG_INFO("sigma", sigma);

            *get<1>(iter) = sigma;
        }
        const typename traitsT::real_type epsilon =
            toml::get<toml::Float>(potent.at("epsilon"));
        MJOLNIR_LOG_INFO("epsilon", epsilon);

        MJOLNIR_LOG_DEBUG("read_global_potential RETURNED");
        return make_unique<ExcludedVolumePotential<traitsT>>(
                epsilon, std::move(params));
    }
    else if(name == "LennardJones")
    {
        toml::Array<toml::Table> param_table =
            toml::get<toml::Array<toml::Table>>(potent.at("parameters"));
        MJOLNIR_LOG_INFO("parameter table size", param_table.size());

        std::vector<typename LennardJonesPotential<traitsT>::parameter_type>
            params(param_table.size());

        for(auto iter = make_zip(param_table.cbegin(), params.begin());
                iter != make_zip(param_table.cend(), params.end()); ++iter)
        {
            const typename traitsT::real_type sigma =
                toml::get<toml::Float>(get<0>(iter)->at("sigma"));
            const typename traitsT::real_type epsilon =
                toml::get<toml::Float>(get<0>(iter)->at("epsilon"));
            MJOLNIR_LOG_INFO("sigma", sigma, "epsilon", epsilon);

            *get<1>(iter) = std::make_pair(sigma, epsilon);
        }

        MJOLNIR_LOG_DEBUG("read_global_potential RETURNED");
        return make_unique<LennardJonesPotential<traitsT>>(
                std::move(params));
    }
    else
        throw std::runtime_error("unknown potential: " + name);
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

    if(name == "Global")
    {
        const typename traitsT::real_type cutoff = pot->max_cutoff_length();
        MJOLNIR_LOG_INFO("potential cutoff length", cutoff);

        typename traitsT::real_type mergin = 1.0;
        try{mergin = toml::get<toml::Float>(potent.at("mergin"));}
        catch(std::exception& except){mergin = 1.0;}
        MJOLNIR_LOG_INFO("mergin rate", mergin);

//         auto space = make_unique<VerletList<traitsT>>(cutoff, cutoff * mergin);
        auto space = make_unique<UnlimitedGridCellList<traitsT>>(
                cutoff, cutoff * mergin);
        MJOLNIR_LOG_DEBUG("space partitioning class created.");

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
        space->set_except(excepts);

        MJOLNIR_LOG_DEBUG("read_global_interaction RETURNED");
        return make_unique<GlobalDistanceInteraction<traitsT>>(std::move(space));
    }
    else
        throw std::runtime_error("unknown interaction: " + name);
}

template<typename traitsT>
GlobalForceField<traitsT>
read_global_force_field(const toml::Array<toml::Table>& gffs)
{
    MJOLNIR_SET_LOGGER("read_toml_file");
    MJOLNIR_LOG_DEBUG("read_global_force_field CALLED");
    GlobalForceField<traitsT> gff;

    for(auto iter = gffs.cbegin(); iter != gffs.cend(); ++iter)
    {
        const std::string potential =
            toml::get<toml::String>(iter->at("potential"));
        MJOLNIR_LOG_INFO("potential name read", potential);

        std::unique_ptr<GlobalPotentialBase<traitsT>> pot =
            read_global_potential<traitsT>(potential, *iter);

        const std::string interaction =
            toml::get<toml::String>(iter->at("interaction"));
        MJOLNIR_LOG_INFO("interaction name read", interaction);

        std::unique_ptr<GlobalInteractionBase<traitsT>> inter =
            read_global_interaction<traitsT>(interaction, pot, *iter);

        gff.emplace(std::move(inter), std::move(pot));
    }
    MJOLNIR_LOG_DEBUG("read_global_force_field RETURNED");
    return gff;
}

} // mjolnir
#endif/* MJOLNIR_IO_TOML_READ_GLOBAL_FORCEFIELD */
