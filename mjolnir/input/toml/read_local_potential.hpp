#ifndef MJOLNIR_IO_TOML_READ_LOCAL_POTENTIAL
#define MJOLNIR_IO_TOML_READ_LOCAL_POTENTIAL
#include <mjolnir/core/LocalForceField.hpp>
#include <mjolnir/potential/HarmonicPotential.hpp>
#include <mjolnir/potential/ClementiDihedralPotential.hpp>
#include <mjolnir/potential/Go1012ContactPotential.hpp>
#include <mjolnir/potential/FlexibleLocalAnglePotential.hpp>
#include <mjolnir/potential/FlexibleLocalDihedralPotential.hpp>
#include <mjolnir/potential/GaussianPotential.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>
#include <toml/toml.hpp>
#include <algorithm>

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

template<typename traitsT>
struct read_local_potential_impl<traitsT, FlexibleLocalAnglePotential<traitsT>>
{
    static FlexibleLocalAnglePotential<traitsT> invoke(const toml::Table& param)
    {
        MJOLNIR_SET_LOGGER("read_toml_file");
        MJOLNIR_LOG_INFO("read flp-angle potential");

        const typename traitsT::real_type k =
            toml::get<toml::Float>(param.at("k"));

        MJOLNIR_LOG_INFO("k", k);

        std::array<typename traitsT::real_type, 10> term1;
        std::array<typename traitsT::real_type, 10> term2;
        const auto t1 = toml::get<toml::Array<toml::Float>>(param.at("term1"));
        const auto t2 = toml::get<toml::Array<toml::Float>>(param.at("term2"));
        if(t1.size() != 10 or t2.size() != 10)
            throw std::invalid_argument("invalid FLP angle term size");

        std::copy(t1.begin(), t1.end(), term1.begin());
        std::copy(t2.begin(), t2.end(), term2.begin());
        MJOLNIR_LOG_DEBUG("read_local_potential RETURNED");

        return FlexibleLocalAnglePotential<traitsT>(k, term1, term2);
    }
};

template<typename traitsT>
struct read_local_potential_impl<traitsT,
                                 FlexibleLocalDihedralPotential<traitsT>>
{
    static FlexibleLocalDihedralPotential<traitsT>
    invoke(const toml::Table& param)
    {
        MJOLNIR_SET_LOGGER("read_toml_file");
        MJOLNIR_LOG_INFO("read flp-dihedral potential");

        const typename traitsT::real_type k =
            toml::get<toml::Float>(param.at("k"));

        MJOLNIR_LOG_INFO("k", k);

        std::array<typename traitsT::real_type, 7> term;
        const auto t = toml::get<toml::Array<toml::Float>>(param.at("term"));
        if(t.size() != 7)
            throw std::invalid_argument("invalid FLP angle term size");

        std::copy(t.begin(), t.end(), term.begin());

        return FlexibleLocalDihedralPotential<traitsT>(k, term);
    }
};

template<typename traitsT>
struct read_local_potential_impl<traitsT,
                                 GaussianPotential<traitsT>>
{
    static GaussianPotential<traitsT>
    invoke(const toml::Table& param)
    {
        MJOLNIR_SET_LOGGER("read_toml_file");
        MJOLNIR_LOG_INFO("read gaussian potential");

        const typename traitsT::real_type e =
            toml::get<toml::Float>(param.at("epsilon"));
        MJOLNIR_LOG_INFO("epsilon", e);

        const typename traitsT::real_type w =
            toml::get<toml::Float>(param.at("w"));

        const typename traitsT::real_type native =
            toml::get<toml::Float>(param.at("native"));
        MJOLNIR_LOG_INFO("native", native);

        return GaussianPotential<traitsT>(e, w, native);
    }
};



} // mjolnir
#endif /* MJOLNIR_IO_TOML_READ_LOCAL_POTENTIAL */
