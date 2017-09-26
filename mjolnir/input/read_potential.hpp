#ifndef MJOLNIR_READ_POTENTIAL
#define MJOLNIR_READ_POTENTIAL
#include <extlib/toml/toml.hpp>
#include <mjolnir/potential/local/HarmonicPotential.hpp>
#include <mjolnir/potential/local/Go1012ContactPotential.hpp>
#include <mjolnir/potential/local/ClementiDihedralPotential.hpp>
#include <mjolnir/potential/local/GaussianPotential.hpp>
#include <mjolnir/potential/local/FlexibleLocalAnglePotential.hpp>
#include <mjolnir/potential/local/FlexibleLocalDihedralPotential.hpp>
#include <mjolnir/potential/global/ExcludedVolumePotential.hpp>
#include <mjolnir/potential/global/LennardJonesPotential.hpp>
#include <mjolnir/potential/global/DebyeHuckelPotential.hpp>
#include <mjolnir/potential/global/ImplicitMembranePotential.hpp>

namespace mjolnir
{

template<typename traitsT, std::size_t N>
std::vector<std::pair<std::array<std::size_t, N>, HarmonicPotential<traitsT>>>
    read_harmonic_potential(const toml::Table& local)
{
    const auto& params = local.at("parameters").cast<toml::value_t::Array>();
    std::vector<std::pair<std::array<std::size_t, N>,
			  HarmonicPotential<traitsT>>> retval; retval.reserve(params.size());

    for(const auto& item : params)
    {
        const auto& parameter = item.cast<toml::value_t::Table>();
        auto indices = toml::get<std::array<std::size_t, N>>(
	    parameter.at("indices"));
        auto r0 = toml::get<typename traitsT::real_type>(parameter.at("native"));
        auto k  = toml::get<typename traitsT::real_type>(parameter.at("k"));

        retval.emplace_back(std::move(indices), HarmonicPotential<traitsT>(k, r0));
    }
    return retval;
}

template<typename traitsT, std::size_t N>
std::vector<std::pair<std::array<std::size_t, N>, Go1012ContactPotential<traitsT>>>
    read_go1012_contact_potential(const toml::Table& local)
{
    const auto& params = local.at("parameters").cast<toml::value_t::Array>();
    std::vector<std::pair<std::array<std::size_t, N>,
			  Go1012ContactPotential<traitsT>>> retval; retval.reserve(params.size());

    for(const auto& item : params)
    {
        const auto& parameter = item.cast<toml::value_t::Table>();
        auto indices = toml::get<std::array<std::size_t, N>>(
	    parameter.at("indices"));
        auto r0 = toml::get<typename traitsT::real_type>(parameter.at("native"));
        auto k  = toml::get<typename traitsT::real_type>(parameter.at("k"));

        retval.emplace_back(std::move(indices),
                            Go1012ContactPotential<traitsT>(k, r0));
    }
    return retval;
}

template<typename traitsT, std::size_t N>
std::vector<std::pair<std::array<std::size_t, N>, GaussianPotential<traitsT>>>
    read_gaussian_potential(const toml::Table& local)
{
    typedef typename traitsT::real_type real_type;
    const auto& params = local.at("parameters").cast<toml::value_t::Array>();
    std::vector<std::pair<std::array<std::size_t, N>,
			  GaussianPotential<traitsT>>> retval; retval.reserve(params.size());

    for(const auto& item : params)
    {
        const auto& parameter = item.cast<toml::value_t::Table>();
        auto indices = toml::get<std::array<std::size_t, N>>(
	    parameter.at("indices"));
        auto native  = toml::get<real_type>(parameter.at("native"));
        auto epsilon = toml::get<real_type>(parameter.at("epsilon"));
        auto w       = toml::get<real_type>(parameter.at("w"));

        retval.emplace_back(std::move(indices),
                            GaussianPotential<traitsT>(epsilon, w, native));
    }
    return retval;
}

template<typename traitsT, std::size_t N>
std::vector<std::pair<std::array<std::size_t, N>,
                      FlexibleLocalAnglePotential<traitsT>>>
    read_flexible_local_angle_potential(const toml::Table& local)
{
    typedef typename traitsT::real_type real_type;

    const auto& params = local.at("parameters").cast<toml::value_t::Array>();
    std::vector<std::pair<std::array<std::size_t, N>,
			  FlexibleLocalAnglePotential<traitsT>>> retval;
    retval.reserve(params.size());

    for(const auto& item : params)
    {
        const auto& parameter = item.cast<toml::value_t::Table>();
        auto indices = toml::get<std::array<std::size_t, N>>(
	    parameter.at("indices"));
        auto k     = toml::get<real_type>(parameter.at("k"));
        auto term1 = toml::get<std::array<real_type, 10>>(parameter.at("term1"));
        auto term2 = toml::get<std::array<real_type, 10>>(parameter.at("term2"));

        retval.emplace_back(std::move(indices),
			    FlexibleLocalAnglePotential<traitsT>(k, term1, term2));
    }
    return retval;
}

template<typename traitsT, std::size_t N>
std::vector<std::pair<std::array<std::size_t, N>,
                      ClementiDihedralPotential<traitsT>>>
    read_clementi_dihedral_potential(const toml::Table& local)
{
    typedef typename traitsT::real_type real_type;
    const auto& params = local.at("parameters").cast<toml::value_t::Array>();
    std::vector<std::pair<std::array<std::size_t, N>,
			  ClementiDihedralPotential<traitsT>>> retval; retval.reserve(params.size());

    for(const auto& item : params)
    {
        const auto& parameter = item.cast<toml::value_t::Table>();
        auto indices = toml::get<std::array<std::size_t, N>>(
	    parameter.at("indices"));
        auto native = toml::get<real_type>(parameter.at("native"));
        auto k1     = toml::get<real_type>(parameter.at("k1"));
        auto k3     = toml::get<real_type>(parameter.at("k3"));

        retval.emplace_back(std::move(indices),
                            ClementiDihedralPotential<traitsT>(k1, k3, native));
    }
    return retval;
}

template<typename traitsT, std::size_t N>
std::vector<std::pair<std::array<std::size_t, N>,
                      FlexibleLocalDihedralPotential<traitsT>>>
    read_flexible_local_dihedral_potential(const toml::Table& local)
{
    typedef typename traitsT::real_type real_type;
    const auto& params = local.at("parameters").cast<toml::value_t::Array>();
    std::vector<std::pair<std::array<std::size_t, N>,
			  FlexibleLocalDihedralPotential<traitsT>>> retval;
    retval.reserve(params.size());

    for(const auto& item : params)
    {
        const auto& parameter = item.cast<toml::value_t::Table>();
        auto indices = toml::get<std::array<std::size_t, N>>(
	    parameter.at("indices"));
        auto k    = toml::get<real_type>(parameter.at("k"));
        auto term = toml::get<std::array<real_type, 7>>(parameter.at("term"));

        retval.emplace_back(std::move(indices),
                            FlexibleLocalDihedralPotential<traitsT>(k, term));
    }
    return retval;
}

template<typename traitsT>
ExcludedVolumePotential<traitsT>
read_excluded_volume_potential(const toml::Table& global)
{
    typedef typename traitsT::real_type real_type;
    real_type eps = toml::get<real_type>(global.at("epsilon"));
    const auto& ps = global.at("parameters").cast<toml::value_t::Array>();

    std::vector<real_type> params; params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = toml::get<std::size_t>(tab.at("index"));
        if(params.size() <= idx)
            params.resize(idx+1, 0.);
        params.at(idx) = toml::get<real_type>(tab.at("sigma"));
    }

    return ExcludedVolumePotential<traitsT>(eps, std::move(params));
}

template<typename traitsT>
LennardJonesPotential<traitsT>
read_lennard_jones_potential(const toml::Table& global)
{
    typedef typename traitsT::real_type real_type;
    const auto& ps = global.at("parameters").cast<toml::value_t::Array>();

    std::vector<std::pair<real_type, real_type>> params; params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = toml::get<std::size_t>(tab.at("index"));
        if(params.size() <= idx)
        {
            const std::pair<real_type, real_type> dummy{0., 0.};
            params.resize(idx+1, dummy);
        }
        params.at(idx) = std::make_pair(toml::get<real_type>(tab.at("sigma")),
                                        toml::get<real_type>(tab.at("epsilon")));
    }

    return LennardJonesPotential<traitsT>(std::move(params));
}

template<typename traitsT>
DebyeHuckelPotential<traitsT>
read_debye_huckel_potential(const toml::Table& global)
{
    typedef typename traitsT::real_type real_type;
    const auto& ps = global.at("parameters").cast<toml::value_t::Array>();

    std::vector<real_type> params; params.reserve(ps.size());
    for(const auto& param : ps)
    {
        const auto& tab = param.cast<toml::value_t::Table>();
        const auto idx = toml::get<std::size_t>(tab.at("index"));
        if(params.size() <= idx)
            params.resize(idx+1, 0.);
        params.at(idx) = toml::get<real_type>(tab.at("charge"));
    }

    return DebyeHuckelPotential<traitsT>(std::move(params));
}

template<typename traitsT>
ImplicitMembranePotential<traitsT>
read_implicit_membrane_potential(const toml::Table& global)
{
    typedef typename traitsT::real_type real_type;
    real_type thickness = toml::get<real_type>(global.at("thickness"));
    real_type interaction_magnitude =
	toml::get<real_type>(global.at("interaction_magnitude"));
    real_type bend = toml::get<real_type>(global.at("bend"));
    const auto& ps = global.at("parameters").cast<toml::value_t::Array>();

    std::vector<real_type> params; params.reserve(ps.size());
    for(const auto& param : ps)
    {
	const auto& tab = param.cast<toml::value_t::Table>();
	const auto idx = toml::get<std::size_t>(tab.at("index"));
	if(params.size() <= idx)
	    params.resize(idx+1, 0.);
	params.at(idx) = toml::get<real_type>(tab.at("hydrophobicity"));
    }
    
    return ImplicitMembranePotential<traitsT>(
	thickness, interaction_magnitude, bend, std::move(params));
}
  
}
#endif // MJOLNIR_READ_POTENTIAL
