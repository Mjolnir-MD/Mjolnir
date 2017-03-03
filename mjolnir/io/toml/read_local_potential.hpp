#ifndef MJOLNIR_IO_TOML_READ_LOCAL_POTENTIAL
#define MJOLNIR_IO_TOML_READ_LOCAL_POTENTIAL
#include <mjolnir/core/LocalForceField.hpp>
#include <mjolnir/potential/HarmonicPotential.hpp>
#include <mjolnir/potential/ClementiDihedralPotential.hpp>
#include <mjolnir/potential/Go1012ContactPotential.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>
#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT, std::size_t N>
struct read_local_indices_impl
{
    static
    std::array<std::size_t, N> invoke(const toml::Table& param)
    {
        const auto indices =
            toml::get<toml::Array<toml::Integer>>(param.at("indices"));

        std::array<std::size_t, N> retval;
        for(std::size_t i=0; i<N; ++i)
            retval[i] = indices.at(i);
        return retval;
    }
};

template<typename traitsT, typename potentialT>
struct read_local_potential_impl;

template<typename traitsT, typename potentialT, std::size_t N>
std::vector<std::pair<std::array<std::size_t, N>, potentialT>>
read_local_potential(const toml::Array<toml::Table>& params)
{
    MJOLNIR_SET_LOGGER("read_toml_file");

    std::vector<std::pair<std::array<std::size_t, N>, potentialT>> retval;
    retval.reserve(params.size());
    for(auto iter = params.cbegin(); iter != params.cend(); ++iter)
    {
        retval.emplace_back(read_local_indices_impl<traitsT, N>::invoke(*iter),
                read_local_potential_impl<traitsT, potentialT>::invoke(*iter));
    }
    return retval;
}

template<typename traitsT>
struct read_local_potential_impl<traitsT, HarmonicPotential<traitsT>>
{
    static HarmonicPotential<traitsT> invoke(const toml::Table& param)
    {
        MJOLNIR_SET_LOGGER("read_toml_file");
        MJOLNIR_LOG_INFO("read harmonic potential");
        const typename traitsT::real_type k =
            toml::get<toml::Float>(param.at("k"));

        MJOLNIR_LOG_INFO("k", k);
        const typename traitsT::real_type r0 =
            toml::get<toml::Float>(param.at("r0"));
        MJOLNIR_LOG_INFO("r0", r0);
        MJOLNIR_LOG_DEBUG("read_local_potential RETURNED");
        return HarmonicPotential<traitsT>(k, r0);
    }
};

template<typename traitsT>
struct read_local_potential_impl<traitsT, ClementiDihedralPotential<traitsT>>
{
    static ClementiDihedralPotential<traitsT> invoke(const toml::Table& param)
    {
        MJOLNIR_SET_LOGGER("read_toml_file");
        MJOLNIR_LOG_INFO("read clementi-dihedral potential");
        const typename traitsT::real_type k1 =
            toml::get<toml::Float>(param.at("k1"));
        MJOLNIR_LOG_INFO("k1", k1);
        const typename traitsT::real_type k3 =
            toml::get<toml::Float>(param.at("k3"));
        MJOLNIR_LOG_INFO("k3", k3);
        const typename traitsT::real_type phi0 =
            toml::get<toml::Float>(param.at("phi0"));
        MJOLNIR_LOG_INFO("phi0", phi0);
        MJOLNIR_LOG_DEBUG("read_local_potential RETURNED");
        return ClementiDihedralPotential<traitsT>(k1, k3, phi0);
    }
};

template<typename traitsT>
struct read_local_potential_impl<traitsT, Go1012ContactPotential<traitsT>>
{
    static Go1012ContactPotential<traitsT> invoke(const toml::Table& param)
    {
        MJOLNIR_SET_LOGGER("read_toml_file");
        MJOLNIR_LOG_INFO("read go-contact potential");
        const typename traitsT::real_type k =
            toml::get<toml::Float>(param.at("k"));
        MJOLNIR_LOG_INFO("k", k);
        const typename traitsT::real_type r0 =
            toml::get<toml::Float>(param.at("r0"));
        MJOLNIR_LOG_INFO("r0", r0);
        MJOLNIR_LOG_DEBUG("read_local_potential RETURNED");
        return Go1012ContactPotential<traitsT>(k, r0);
    }
};


} // mjolnir
#endif /* MJOLNIR_IO_TOML_READ_LOCAL_POTENTIAL */
