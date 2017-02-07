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

template<typename traitsT>
std::unique_ptr<LocalPotentialBase<traitsT>>
read_harmonic_potential(const toml::Table& param)
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
    return make_unique<HarmonicPotential<traitsT>>(k, r0);
}

template<typename traitsT>
std::unique_ptr<LocalPotentialBase<traitsT>>
read_clementi_dihedral_potential(const toml::Table& param)
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
    return make_unique<ClementiDihedralPotential<traitsT>>(k1, k3, phi0);
}

template<typename traitsT>
std::unique_ptr<LocalPotentialBase<traitsT>>
read_go_1012_contact_potential(const toml::Table& param)
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
    return make_unique<Go1012ContactPotential<traitsT>>(k, r0);
}

template<typename traitsT>
std::unique_ptr<LocalPotentialBase<traitsT>>
read_local_potential(const std::string& name, const toml::Table& param)
{
    MJOLNIR_SET_LOGGER("read_toml_file");
    MJOLNIR_LOG_DEBUG("read_local_potential CALLED");
    MJOLNIR_LOG_INFO("local potential name", name);

    if(name == "Harmonic")
        return read_harmonic_potential<traitsT>(param);
    else if(name == "ClementiDihedral")
        return read_clementi_dihedral_potential<traitsT>(param);
    else if(name == "Go1012Contact")
        return read_go_1012_contact_potential<traitsT>(param);
    else
        throw std::runtime_error("unknown potential: " + name);
}

template<typename traitsT, std::size_t N>
typename LocalForceField<traitsT>::template potential_array<N>
read_local_potential_array(
        const std::string& potential, const toml::Array<toml::Table>& params)
{
    MJOLNIR_SET_LOGGER("read_toml_file");
    typename LocalForceField<traitsT>::template potential_array<N> pots;

    for(auto p = params.cbegin(); p != params.cend(); ++p)
    {
        std::unique_ptr<LocalPotentialBase<traitsT>> pot =
            read_local_potential<traitsT>(potential, *p);
        const toml::Array<toml::Integer> indices =
            toml::get<toml::Array<toml::Integer>>(p->at("indices"));
        MJOLNIR_LOG_INFO("indices.size()", indices.size());
        if(indices.size() != N)
            throw std::runtime_error("bond index must be specified");

        std::array<std::size_t, N> idx;
        for(std::size_t i=0; i < N; ++i)
        {
            idx[i] = indices[i];
            MJOLNIR_LOG_DEBUG("index at ", i, "is", indices.at(i));
        }
        pots.emplace_back(std::move(idx), std::move(pot));
    }
    return pots;
}

} // mjolnir
#endif /* MJOLNIR_IO_TOML_READ_LOCAL_POTENTIAL */
