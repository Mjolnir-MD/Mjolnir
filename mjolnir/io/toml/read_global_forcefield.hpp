#ifndef MJOLNIR_IO_TOML_READ_GLOBAL_FORCEFIELD
#define MJOLNIR_IO_TOML_READ_GLOBAL_FORCEFIELD
#include <mjolnir/core/GlobalDistanceInteraction.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <mjolnir/potential/ExcludedVolumePotential.hpp>
#include <mjolnir/potential/LennardJonesPotential.hpp>
#include <mjolnir/util/zip_iterator.hpp>
#include <mjolnir/util/make_zip.hpp>
#include <mjolnir/util/make_unique.hpp>

#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<GlobalPotentialBase<traitsT>>
read_global_potential(const std::string& name, const toml::Table& potent)
{
    if(name == "ExcludedVolume")
    {
        toml::Array<toml::Table> param_table =
            toml::get<toml::Array<toml::Table>>(potent.at("parameters"));
        std::vector<typename ExcludedVolumePotential<traitsT>::parameter_type>
            params(param_table.size());
        for(auto iter = make_zip(param_table.cbegin(), params.begin());
                iter != make_zip(param_table.cend(), params.end()); ++iter)
        {
            const typename traitsT::real_type sigma = 
                toml::get<toml::Float>(get<0>(iter)->at("sigma"));

            *get<1>(iter) = sigma;
        }
        const typename traitsT::real_type epsilon = 
            toml::get<toml::Float>(potent.at("epsilon"));

        return make_unique<ExcludedVolumePotential<traitsT>>(
                epsilon, std::move(params));
    }
    else if(name == "LennardJones")
    {
        toml::Array<toml::Table> param_table =
            toml::get<toml::Array<toml::Table>>(potent.at("parameters"));
        std::vector<typename LennardJonesPotential<traitsT>::parameter_type>
            params(param_table.size());

        for(auto iter = make_zip(param_table.cbegin(), params.begin());
                iter != make_zip(param_table.cend(), params.end()); ++iter)
        {
            const typename traitsT::real_type sigma = 
                toml::get<toml::Float>(get<0>(iter)->at("sigma"));
            const typename traitsT::real_type epsilon = 
                toml::get<toml::Float>(get<0>(iter)->at("epsilon"));

            *get<1>(iter) = std::make_pair(sigma, epsilon);
        }

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
    if(name == "Global")
    {
        const typename traitsT::real_type cutoff = pot->max_cutoff_length();

        typename traitsT::real_type mergin = 1.0;
        try{mergin = toml::get<toml::Float>(potent.at("mergin"));}
        catch(toml::exception& except){mergin = 1.0;}

        auto space = make_unique<VerletList<traitsT>>(cutoff, cutoff * mergin);

        std::vector<std::vector<std::size_t>> excepts;
        auto elist = toml::get<toml::Array<toml::Array<toml::Integer>>>(
                potent.at("excepts"));
        for(auto iter = elist.cbegin(); iter != elist.cend(); ++iter)
        {
            std::vector<std::size_t> l; l.reserve(iter->size());
            for(auto j = iter->cbegin(); j != iter->cend(); ++j)
                l.push_back(static_cast<std::size_t>(*j));
            excepts.emplace_back(l);
        }
        space->set_except(excepts);

        return make_unique<GlobalDistanceInteraction<traitsT>>(std::move(space));
    }
    else
        throw std::runtime_error("unknown interaction: " + name);
}

template<typename traitsT>
GlobalForceField<traitsT>
read_global_force_field(const toml::Array<toml::Table>& gffs)
{
    GlobalForceField<traitsT> gff;

    for(auto iter = gffs.cbegin(); iter != gffs.cend(); ++iter)
    {
        const std::string potential = 
            toml::get<toml::String>(iter->at("potential"));
        std::unique_ptr<GlobalPotentialBase<traitsT>> pot = 
            read_global_potential<traitsT>(potential, *iter);

        const std::string interaction = 
            toml::get<toml::String>(iter->at("interaction"));
        std::unique_ptr<GlobalInteractionBase<traitsT>> inter = 
            read_global_interaction<traitsT>(interaction, pot, *iter);

        gff.emplace(std::move(inter), std::move(pot));
    }
    return gff;
}

} // mjolnir
#endif/* MJOLNIR_IO_TOML_READ_GLOBAL_FORCEFIELD */
