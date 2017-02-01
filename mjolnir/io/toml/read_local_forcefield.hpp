#ifndef MJOLNIR_IO_TOML_READ_LOCAL_FORCEFIELD
#define MJOLNIR_IO_TOML_READ_LOCAL_FORCEFIELD
#include <mjolnir/core/BondLengthInteraction.hpp>
#include <mjolnir/core/BondAngleInteraction.hpp>
#include <mjolnir/core/DihedralAngleInteraction.hpp>
#include <mjolnir/core/LocalForceField.hpp>
#include <mjolnir/potential/HarmonicPotential.hpp>
#include <mjolnir/potential/ClementiDihedralPotential.hpp>
#include <mjolnir/potential/Go1012ContactPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

#include <toml/toml.hpp>

namespace mjolnir
{

template<typename traitsT>
std::unique_ptr<LocalPotentialBase<traitsT>>
read_local_potential(const std::string& name, const toml::Table& param)
{
    MJOLNIR_SET_LOGGER("read_toml_file");
    MJOLNIR_LOG_DEBUG("read_local_potential CALLED");
    MJOLNIR_LOG_INFO("local potential name", name);

    if(name == "Harmonic")
    {
        const typename traitsT::real_type k =
            toml::get<toml::Float>(param.at("k"));
        MJOLNIR_LOG_INFO("k", k);
        const typename traitsT::real_type r0 =
            toml::get<toml::Float>(param.at("r0"));
        MJOLNIR_LOG_INFO("r0", r0);
        MJOLNIR_LOG_DEBUG("read_local_potential RETURNED");
        return make_unique<HarmonicPotential<traitsT>>(k, r0);
    }
    else if(name == "ClementiDihedral")
    {
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
    else if(name == "Go1012Contact")
    {
        const typename traitsT::real_type k =
            toml::get<toml::Float>(param.at("k"));
        MJOLNIR_LOG_INFO("k", k);
        const typename traitsT::real_type r0 =
            toml::get<toml::Float>(param.at("r0"));
        MJOLNIR_LOG_INFO("r0", r0);
        MJOLNIR_LOG_DEBUG("read_local_potential RETURNED");
        return make_unique<Go1012ContactPotential<traitsT>>(k, r0);
    }
    else
        throw std::runtime_error("unknown potential: " + name);
}

template<typename traitsT>
LocalForceField<traitsT>
read_local_force_field(const toml::Array<toml::Table>& lffs)
{
    MJOLNIR_SET_LOGGER("read_toml_file");
    MJOLNIR_LOG_DEBUG("read_local_force_field CALLED");
    LocalForceField<traitsT> lff;
    for(auto iter = lffs.cbegin(); iter != lffs.cend(); ++iter)
    {
        const std::string potential =
            toml::get<toml::String>(iter->at("potential"));
        MJOLNIR_LOG_INFO("potential name", potential);

        const toml::Array<toml::Table> params =
            toml::get<toml::Array<toml::Table>>(iter->at("parameters"));
        MJOLNIR_LOG_INFO("parameter table size", params.size());

        const std::string interaction =
            toml::get<toml::String>(iter->at("interaction"));
        MJOLNIR_LOG_INFO("interaction name", interaction);

        if(interaction == "BondLength")
        {
            for(auto p = params.cbegin(); p != params.cend(); ++p)
            {
                std::unique_ptr<LocalPotentialBase<traitsT>> pot =
                    read_local_potential<traitsT>(potential, *p);
                const toml::Array<toml::Integer> indices =
                    toml::get<toml::Array<toml::Integer>>(p->at("indices"));
                MJOLNIR_LOG_INFO("indices.size()", indices.size());
                if(indices.size() != 2)
                    throw std::runtime_error("bond index must be specified");

                MJOLNIR_LOG_INFO("index", indices.at(0), indices.at(1));
                lff.emplace_bond(indices.at(0), indices.at(1), std::move(pot));
            }
        }
        else if(interaction == "BondAngle")
        {
            for(auto p = params.cbegin(); p != params.cend(); ++p)
            {
                std::unique_ptr<LocalPotentialBase<traitsT>> pot =
                    read_local_potential<traitsT>(potential, *p);
                const toml::Array<toml::Integer> indices =
                    toml::get<toml::Array<toml::Integer>>(p->at("indices"));
                MJOLNIR_LOG_INFO("indices.size()", indices.size());
                if(indices.size() != 3)
                    throw std::runtime_error("angle index must be specified");
                MJOLNIR_LOG_INFO("index", indices.at(0), indices.at(1),
                                 indices.at(2));
                lff.emplace_angle(indices.at(0), indices.at(1), indices.at(2),
                                  std::move(pot));
            }
        }
        else if(interaction == "DihedralAngle")
        {
            for(auto p = params.cbegin(); p != params.cend(); ++p)
            {
                std::unique_ptr<LocalPotentialBase<traitsT>> pot =
                    read_local_potential<traitsT>(potential, *p);
                const toml::Array<toml::Integer> indices =
                    toml::get<toml::Array<toml::Integer>>(p->at("indices"));
                MJOLNIR_LOG_INFO("indices.size()", indices.size());
                if(indices.size() != 4)
                    throw std::runtime_error("dihedral index must be specified");
                MJOLNIR_LOG_INFO("index", indices.at(0), indices.at(1),
                                 indices.at(2), indices.at(3));
                lff.emplace_dihedral(indices.at(0), indices.at(1), indices.at(2),
                        indices.at(3), std::move(pot));
            }
        }
        else
            throw std::runtime_error("unknown interaction: " + interaction);
    }
    MJOLNIR_LOG_DEBUG("read_local_force_field RETUENED");
    return lff;
}


} // mjolnir
#endif /*MJOLNIR_IO_TOML_READ_LOCAL_FORCEFIELD */
