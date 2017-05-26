#ifndef MJOLNIR_IO_TOML_READ_GLOBAL_POTENTIAL
#define MJOLNIR_IO_TOML_READ_GLOBAL_POTENTIAL
#include <mjolnir/potential/ExcludedVolumePotential.hpp>
#include <mjolnir/potential/LennardJonesPotential.hpp>
#include <mjolnir/util/zip_iterator.hpp>
#include <mjolnir/util/make_zip.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>
#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT, typename potentialT>
struct read_global_potential_impl;

template<typename traitsT, typename potentialT>
potentialT
read_global_potential(const toml::Table& pot)
{
    return read_global_potential_impl<traitsT, potentialT>::invoke(pot);
}

template<typename traitsT>
struct read_global_potential_impl<traitsT, ExcludedVolumePotential<traitsT>>
{
    static ExcludedVolumePotential<traitsT>
    invoke(const toml::Table& potent)
    {
        MJOLNIR_SET_LOGGER("read_toml_file");
        MJOLNIR_LOG_DEBUG("read_excluded_volume CALLED");

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
        return ExcludedVolumePotential<traitsT>(
                std::move(epsilon), std::move(params));
    }
};

template<typename traitsT>
struct read_global_potential_impl<traitsT, LennardJonesPotential<traitsT>>
{
    static LennardJonesPotential<traitsT>
    invoke(const toml::Table& potent)
    {
        MJOLNIR_SET_LOGGER("read_toml_file");
        MJOLNIR_LOG_DEBUG("read_lennard_jones CALLED");
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
        return LennardJonesPotential<traitsT>(std::move(params));
    }
};

} // mjolnir
#endif /* MJOLNIR_IO_TOML_READ_GLOBAL_POTENTIAL */
