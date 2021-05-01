#ifndef MJOLNIR_INPUT_READ_LOCAL_POTENTIAL_HPP
#define MJOLNIR_INPUT_READ_LOCAL_POTENTIAL_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/input/utility.hpp>

#include <mjolnir/forcefield/local/HarmonicPotential.hpp>
#include <mjolnir/forcefield/local/GoContactPotential.hpp>
#include <mjolnir/forcefield/local/GoContactAttractivePotential.hpp>
#include <mjolnir/forcefield/local/GoContactRepulsivePotential.hpp>
#include <mjolnir/forcefield/local/ClementiDihedralPotential.hpp>
#include <mjolnir/forcefield/local/GaussianPotential.hpp>
#include <mjolnir/forcefield/local/PeriodicGaussianPotential.hpp>
#include <mjolnir/forcefield/local/CosinePotential.hpp>
#include <mjolnir/forcefield/local/SumLocalPotential.hpp>
#include <mjolnir/forcefield/local/UniformPotential.hpp>
#include <mjolnir/forcefield/local/WormLikeChainPotential.hpp>
#include <mjolnir/forcefield/local/WormLikeChainOffsetPotential.hpp>
#include <mjolnir/forcefield/MultipleBasin/MBasinAttractivePotential.hpp>
#include <mjolnir/forcefield/MultipleBasin/MBasinRepulsivePotential.hpp>
#include <mjolnir/forcefield/FLP/FlexibleLocalAnglePotential.hpp>
#include <mjolnir/forcefield/FLP/FlexibleLocalDihedralPotential.hpp>
#include <mjolnir/forcefield/3SPN2/ThreeSPN2BondPotential.hpp>

#include <mjolnir/core/Topology.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// this reads...
// parameters = [
// {indices = [1, 2], k = 10.0, v0 = 1.0} <- a portion of this table, k and v0.
// ]
template<typename realT>
HarmonicPotential<realT>
read_harmonic_potential(const toml::value& param, const toml::value& env)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    check_keys_available(param, {"indices", "k", "v0", "offset"});

    const auto k  = find_parameter<real_type>(param, env, "k" );
    const auto v0 = find_parameter<real_type>(param, env, "v0");

    MJOLNIR_LOG_INFO("HarmonicPotential = {v0 = ", v0, ", k = ", k, '}');
    return HarmonicPotential<realT>(k, v0);
}

template<typename realT>
GoContactPotential<realT>
read_go_contact_potential(const toml::value& param, const toml::value& env)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    check_keys_available(param, {"indices", "k", "v0", "offset"});

    const auto k  = find_parameter<real_type>(param, env, "k");
    const auto v0 = find_parameter<real_type>(param, env, "v0");

    MJOLNIR_LOG_INFO("GoContactPotential = {v0 = ", v0, ", k = ", k, '}');
    return GoContactPotential<realT>(k, v0);
}

template<typename realT>
GoContactAttractivePotential<realT>
read_go_contact_attractive_potential(
        const toml::value& param, const toml::value& env)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    check_keys_available(param, {"indices", "k", "v0", "offset"});

    const auto k  = find_parameter<real_type>(param, env, "k");
    const auto v0 = find_parameter<real_type>(param, env, "v0");

    MJOLNIR_LOG_INFO("GoContactAttractivePotential = {v0 = ", v0, ", k = ", k, '}');
    return GoContactAttractivePotential<realT>(k, v0);
}

template<typename realT>
GoContactRepulsivePotential<realT>
read_go_contact_repulsive_potential(
        const toml::value& param, const toml::value& env)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    check_keys_available(param, {"indices", "k", "v0", "offset"});

    const auto k  = find_parameter<real_type>(param, env, "k");
    const auto v0 = find_parameter<real_type>(param, env, "v0");

    MJOLNIR_LOG_INFO("GoContactRepulsivePotential = {v0 = ", v0, ", k = ", k, '}');
    return GoContactRepulsivePotential<realT>(k, v0);
}

template<typename realT>
MBasinAttractivePotential<realT>
read_mbasin_attractive_potential(
        const toml::value& param, const toml::value& env)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    check_keys_available(param, {"indices", "k", "v0", "offset"});

    const auto k  = find_parameter<real_type>(param, env, "k");
    const auto v0 = find_parameter<real_type>(param, env, "v0");

    MJOLNIR_LOG_INFO("MBasinAttractivePotential = {v0 = ", v0, ", k = ", k, '}');
    return MBasinAttractivePotential<realT>(k, v0);
}

template<typename realT>
MBasinRepulsivePotential<realT>
read_mbasin_repulsive_potential(
        const toml::value& param, const toml::value& env)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    check_keys_available(param, {"indices", "k", "v0", "offset"});

    const auto k  = find_parameter<real_type>(param, env, "k");
    const auto v0 = find_parameter<real_type>(param, env, "v0");

    MJOLNIR_LOG_INFO("MBasinRepulsivePotential = {v0 = ", v0, ", k = ", k, '}');
    return MBasinRepulsivePotential<realT>(k, v0);
}


template<typename realT>
GaussianPotential<realT>
read_gaussian_potential(const toml::value& param, const toml::value& env)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    check_keys_available(param, {"indices", "k", "v0", "sigma", u8"σ", "offset"});

    const auto v0    = find_parameter<real_type>(param, env, "v0");
    const auto k     = find_parameter<real_type>(param, env, "k");
    const auto sigma = find_parameter<real_type>(param, env, "sigma", u8"σ");

    MJOLNIR_LOG_INFO("GaussianPotential = {v0 = ", v0, ", k = ", k,
                     ", sigma = ", sigma, '}');
    return GaussianPotential<realT>(k, sigma, v0);
}

template<typename realT>
PeriodicGaussianPotential<realT>
read_periodic_gaussian_potential(const toml::value& param, const toml::value& env)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    check_keys_available(param, {"indices", "k", "v0", "sigma", u8"σ", "offset"});

    const auto v0    = find_parameter<real_type>(param, env, "v0");
    const auto k     = find_parameter<real_type>(param, env, "k");
    const auto sigma = find_parameter<real_type>(param, env, "sigma", u8"σ");

    MJOLNIR_LOG_INFO("PeriodicGaussianPotential = {v0 = ", v0, ", k = ", k,
                     ", sigma = ", sigma, '}');
    return PeriodicGaussianPotential<realT>(k, sigma, v0);
}

template<typename realT>
FlexibleLocalAnglePotential<realT>
read_flexible_local_angle_potential(const toml::value& param, const toml::value& env)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    check_keys_available(param, {"indices", "k", "y", "d2y", "x", "offset"});
    const auto k     = find_parameter<real_type                >(param, env, "k");
    const auto term1 = find_parameter<std::array<real_type, 10>>(param, env, "y");
    const auto term2 = find_parameter<std::array<real_type, 10>>(param, env, "d2y");

    if(param.as_table().count("x") == 1)
    {
        const auto xs =
            find_parameter<std::array<real_type, 10>>(param, env, "x");
        MJOLNIR_LOG_INFO("FlexibleLocalAngle = {k = ", k, ", x = ", xs,
                         ", y = ", term1, ", d2y = ", term2, '}');
        try
        {
            return FlexibleLocalAnglePotential<realT>(k, xs, term1, term2);
        }
        catch (const std::runtime_error&)
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_flexible_local_angle_potential: data points "
                "are not evenly distributed", find_parameter<toml::value>(
                param, env, "x"), "here"));
        }
    }
    else
    {
        MJOLNIR_LOG_INFO("FlexibleLocalAngle = {k = ", k,
                         ", y = ", term1, ", d2y = ", term2, '}');
        return FlexibleLocalAnglePotential<realT>(k, term1, term2);
    }
}

template<typename realT>
ClementiDihedralPotential<realT>
read_clementi_dihedral_potential(const toml::value& param, const toml::value& env)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    check_keys_available(param, {"indices", "k1", "k3", "v0", "offset"});

    const auto v0 = find_parameter<real_type>(param, env, "v0");
    const auto k1 = find_parameter<real_type>(param, env, "k1");
    const auto k3 = find_parameter<real_type>(param, env, "k3");

    MJOLNIR_LOG_INFO("ClementiDihedral = {v0 = ", v0,
                     ", k1 = ", k1, ", k3 = ", k3, '}');
    return ClementiDihedralPotential<realT>(k1, k3, v0);
}

template<typename realT>
FlexibleLocalDihedralPotential<realT>
read_flexible_local_dihedral_potential(const toml::value& param, const toml::value& env)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    check_keys_available(param, {"indices", "k", "coef", "offset"});

    auto k    = find_parameter<real_type               >(param, env, "k");
    auto coef = find_parameter<std::array<real_type, 7>>(param, env, "coef");

    MJOLNIR_LOG_INFO("FlexibleLocalDihedral = {k = ", k, ", coef = ", coef, '}');
    return FlexibleLocalDihedralPotential<realT>(k, coef);
}

template<typename realT>
CosinePotential<realT>
read_cosine_potential(const toml::value& param, const toml::value& env)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    check_keys_available(param, {"indices", "k", "n", "v0", "offset"});

    auto k  = find_parameter<real_type   >(param, env, "k");
    auto n  = find_parameter<std::int32_t>(param, env, "n");
    auto v0 = find_parameter<real_type   >(param, env, "v0");

    MJOLNIR_LOG_INFO("CosinePotential = {k = ", k, ", n = ", n, ", v0 = ", v0, '}');
    return CosinePotential<realT>(k, n, v0);
}

template<typename realT>
UniformPotential<realT>
read_uniform_potential(const toml::value& param, const toml::value& env)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    check_keys_available(param, {"indices", "k", "offset"});

    auto k = find_parameter<real_type>(param, env, "k");

    MJOLNIR_LOG_INFO("UniformPotential = {k = ", k, '}');
    return UniformPotential<realT>(k);
}

template<typename realT>
WormLikeChainPotential<realT>
read_worm_like_chain_potential(const toml::value& param, const toml::value& env)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    check_keys_available(param, {"indices", "p", "lc", "offset"});

    auto p  = find_parameter<real_type>(param, env, "p");
    auto lc = find_parameter<real_type>(param, env, "lc");

    MJOLNIR_LOG_INFO("WormLikeChainPotential = {p = ", p, ", lc = ", lc,"}");
    return WormLikeChainPotential<realT>(p, lc);
}

template<typename realT>
WormLikeChainOffsetPotential<realT>
read_worm_like_chain_offset_potential(const toml::value& param, const toml::value& env)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    check_keys_available(param, {"indices", "p", "lc", "offset", "l0"});

    auto p  = find_parameter<real_type>(param, env, "p");
    auto lc = find_parameter<real_type>(param, env, "lc");
    auto l0 = find_parameter<real_type>(param, env, "l0");

    MJOLNIR_LOG_INFO("WormLikeChainOffsetPotential = {p = ", p, ", lc = ", lc, " l0 = ", l0, "}");
    return WormLikeChainOffsetPotential<realT>(p, lc, l0);
}


template<typename realT>
ThreeSPN2BondPotential<realT>
read_3spn2_bond_potential(const toml::value& param, const toml::value& env)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    check_keys_available(param, {"indices", "k", "v0", "offset"});

    auto k  = find_parameter<real_type>(param, env, "k");
    auto v0 = find_parameter<real_type>(param, env, "v0");

    MJOLNIR_LOG_INFO("ThreeSPN2BondPotential = {k = ", k, ", v0 = ", v0, '}');
    return ThreeSPN2BondPotential<realT>(k, v0);
}

// ----------------------------------------------------------------------------
// utility function to read local potentials
// ----------------------------------------------------------------------------

// potential_class -> reading_function adapter to implement a function
// `read_local_potential` that reads indices and call read_*_potential.
// local potential requires not only the function parameter but also the indices
// of particles. to split them, mjolnir has both `read_local_potential` and
// `read_(harmonic|gaussian|...)_potential`.
namespace detail
{
template<typename potentialT> struct read_local_potential_impl;

template<typename realT>
struct read_local_potential_impl<HarmonicPotential<realT>>
{
    static HarmonicPotential<realT>
    invoke(const toml::value& param, const toml::value& env)
    {
        return read_harmonic_potential<realT>(param, env);
    }
};
template<typename realT>
struct read_local_potential_impl<GoContactPotential<realT>>
{
    static GoContactPotential<realT>
    invoke(const toml::value& param, const toml::value& env)
    {
        return read_go_contact_potential<realT>(param, env);
    }
};
template<typename realT>
struct read_local_potential_impl<GoContactAttractivePotential<realT>>
{
    static GoContactAttractivePotential<realT>
    invoke(const toml::value& param, const toml::value& env)
    {
        return read_go_contact_attractive_potential<realT>(param, env);
    }
};
template<typename realT>
struct read_local_potential_impl<GoContactRepulsivePotential<realT>>
{
    static GoContactRepulsivePotential<realT>
    invoke(const toml::value& param, const toml::value& env)
    {
        return read_go_contact_repulsive_potential<realT>(param, env);
    }
};
template<typename realT>
struct read_local_potential_impl<MBasinRepulsivePotential<realT>>
{
    static MBasinRepulsivePotential<realT>
    invoke(const toml::value& param, const toml::value& env)
    {
        return read_mbasin_repulsive_potential<realT>(param, env);
    }
};
template<typename realT>
struct read_local_potential_impl<MBasinAttractivePotential<realT>>
{
    static MBasinAttractivePotential<realT>
    invoke(const toml::value& param, const toml::value& env)
    {
        return read_mbasin_attractive_potential<realT>(param, env);
    }
};
template<typename realT>
struct read_local_potential_impl<GaussianPotential<realT>>
{
    static GaussianPotential<realT>
    invoke(const toml::value& param, const toml::value& env)
    {
        return read_gaussian_potential<realT>(param, env);
    }
};
template<typename realT>
struct read_local_potential_impl<PeriodicGaussianPotential<realT>>
{
    static PeriodicGaussianPotential<realT>
    invoke(const toml::value& param, const toml::value& env)
    {
        return read_periodic_gaussian_potential<realT>(param, env);
    }
};
template<typename realT>
struct read_local_potential_impl<FlexibleLocalAnglePotential<realT>>
{
    static FlexibleLocalAnglePotential<realT>
    invoke(const toml::value& param, const toml::value& env)
    {
        return read_flexible_local_angle_potential<realT>(param, env);
    }
};
template<typename realT>
struct read_local_potential_impl<ClementiDihedralPotential<realT>>
{
    static ClementiDihedralPotential<realT>
    invoke(const toml::value& param, const toml::value& env)
    {
        return read_clementi_dihedral_potential<realT>(param, env);
    }
};
template<typename realT>
struct read_local_potential_impl<FlexibleLocalDihedralPotential<realT>>
{
    static FlexibleLocalDihedralPotential<realT>
    invoke(const toml::value& param, const toml::value& env)
    {
        return read_flexible_local_dihedral_potential<realT>(param, env);
    }
};
template<typename realT>
struct read_local_potential_impl<CosinePotential<realT>>
{
    static CosinePotential<realT>
    invoke(const toml::value& param, const toml::value& env)
    {
        return read_cosine_potential<realT>(param, env);
    }
};
template<typename realT>
struct read_local_potential_impl<UniformPotential<realT>>
{
    static UniformPotential<realT>
    invoke(const toml::value& param, const toml::value& env)
    {
        return read_uniform_potential<realT>(param, env);
    }
};
template<typename realT>
struct read_local_potential_impl<WormLikeChainPotential<realT>>
{
    static WormLikeChainPotential<realT>
    invoke(const toml::value& param, const toml::value& env)
    {
        return read_worm_like_chain_potential<realT>(param, env);
    }
};
template<typename realT>
struct read_local_potential_impl<WormLikeChainOffsetPotential<realT>>
{
    static WormLikeChainOffsetPotential<realT>
    invoke(const toml::value& param, const toml::value& env)
    {
        return read_worm_like_chain_offset_potential<realT>(param, env);
    }
};
template<typename realT>
struct read_local_potential_impl<ThreeSPN2BondPotential<realT>>
{
    static ThreeSPN2BondPotential<realT>
    invoke(const toml::value& param, const toml::value& env)
    {
        return read_3spn2_bond_potential<realT>(param, env);
    }
};
} // namespace detail

// this function reads particle indices on which the potential will be applied
// and returns pairs of [indices, potential parameters].
template<std::size_t N, typename potentialT>
std::vector<std::pair<std::array<std::size_t, N>, potentialT>>
read_local_potential(const toml::value& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    MJOLNIR_LOG_INFO("as ", N, "-body interaction");

    using indices_t                = std::array<std::size_t, N>;
    using indices_potential_pair_t = std::pair<indices_t, potentialT>;

    const auto& params = toml::find<toml::array>(local, "parameters");
    MJOLNIR_LOG_NOTICE("-- ", params.size(), " interactions are found.");

    const auto& env = local.contains("env") ? local.at("env") : toml::value{};

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());
    for(const auto& item : params)
    {
        const auto offset = find_parameter_or<toml::value>(item, env, "offset", toml::value(0));

        auto indices = find_parameter<indices_t>(item, env, "indices");

        if(offset.is_array())
        {
            if(offset.size() != indices.size())
            {
                throw_exception<std::runtime_error>("[error] "
                    "mjolnir::read_local_potentials:",
                    "invalid size of offset array", offset.size(), "here",
                    "the size of offset array must much to the size of indices array.",
                    "expected value is ", indices.size(), ".");
            }

            for(std::size_t i=0; i<indices.size(); ++i)
            {
                indices[i] += offset[i].as_integer();
            }
        }
        else
        {
            for(auto& i : indices) {i += offset.as_integer();}
        }

        MJOLNIR_LOG_INFO_NO_LF("idxs = ", indices, ", ");

        retval.emplace_back(indices,
            detail::read_local_potential_impl<potentialT>::invoke(item, env));
    }
    return retval;
}

template<std::size_t N, typename realT,
         template <typename Real> class potential1T,
         template <typename Real> class potential2T
         >
std::vector<std::pair<std::array<std::size_t, N>,
                      SumLocalPotential<realT, potential1T, potential2T>>>
read_local_potentials(const toml::value& local,
        const std::string& pot1_name, const std::string& pot2_name)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    MJOLNIR_LOG_INFO("as ", N, "-body interaction");

    using indices_type   = std::array<std::size_t, N>;
    using potential_type = SumLocalPotential<realT, potential1T, potential2T>;
    using indices_potential_pair_type = std::pair<indices_type, potential_type>;

    const auto& params = toml::find<toml::array>(local, "parameters");
    MJOLNIR_LOG_NOTICE("-- ", params.size(), " interactions are found.");

    const auto& env = local.contains("env") ? local.at("env") : toml::value{};

    std::vector<indices_potential_pair_type> retval;
    retval.reserve(params.size());
    for(const auto& item : params)
    {
        using real_type        = realT;
        using potential_1_type = potential1T<real_type>;
        using potential_2_type = potential2T<real_type>;

        const auto offset = find_parameter_or<toml::value>(item, env, "offset", toml::value(0));

        auto indices = find_parameter<indices_type>(item, env, "indices");

        if(offset.is_array())
        {
            if(offset.size() != indices.size())
            {
                throw_exception<std::runtime_error>("[error] "
                    "mjolnir::read_local_potentials:",
                    "invalid size of offset array", offset.size(), "here",
                    "the size of offset array must much to the size of indices array.",
                    "expected value is ", indices.size(), ".");
            }

            for(std::size_t i=0; i<indices.size(); ++i)
            {
                indices[i] += offset[i].as_integer();
            }
        }
        else
        {
            for(auto& i : indices) {i += offset.as_integer();}
        }


        MJOLNIR_LOG_INFO("idxs = ", indices);

        const auto& pot1 = find_parameter<toml::value>(item, env, pot1_name);
        const auto& pot2 = find_parameter<toml::value>(item, env, pot2_name);

        retval.emplace_back(indices, potential_type(
            detail::read_local_potential_impl<potential_1_type>::invoke(pot1, env),
            detail::read_local_potential_impl<potential_2_type>::invoke(pot2, env))
        );
    }
    return retval;
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template std::vector<std::pair<std::array<std::size_t, 2>, HarmonicPotential<double>>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 2>, HarmonicPotential<float >>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 3>, HarmonicPotential<double>>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 3>, HarmonicPotential<float >>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 4>, HarmonicPotential<double>>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 4>, HarmonicPotential<float >>> read_local_potential(const toml::value& local);

extern template std::vector<std::pair<std::array<std::size_t, 2>, GoContactPotential<double>>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 2>, GoContactPotential<float >>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 2>, GoContactAttractivePotential<double>>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 2>, GoContactAttractivePotential<float >>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 2>, GoContactRepulsivePotential<double>>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 2>, GoContactRepulsivePotential<float >>> read_local_potential(const toml::value& local);

extern template std::vector<std::pair<std::array<std::size_t, 4>, ClementiDihedralPotential<double>>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 4>, ClementiDihedralPotential<float >>> read_local_potential(const toml::value& local);

extern template std::vector<std::pair<std::array<std::size_t, 2>, GaussianPotential<double>>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 2>, GaussianPotential<float >>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 3>, GaussianPotential<double>>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 3>, GaussianPotential<float >>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 4>, PeriodicGaussianPotential<double>>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 4>, PeriodicGaussianPotential<float >>> read_local_potential(const toml::value& local);

extern template std::vector<std::pair<std::array<std::size_t, 3>, FlexibleLocalAnglePotential<double>>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 3>, FlexibleLocalAnglePotential<float >>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 4>, FlexibleLocalDihedralPotential<double>>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 4>, FlexibleLocalDihedralPotential<float >>> read_local_potential(const toml::value& local);

extern template std::vector<std::pair<std::array<std::size_t, 2>, ThreeSPN2BondPotential<double>>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 2>, ThreeSPN2BondPotential<float >>> read_local_potential(const toml::value& local);

extern template std::vector<std::pair<std::array<std::size_t, 4>, CosinePotential<double>>> read_local_potential(const toml::value& local);
extern template std::vector<std::pair<std::array<std::size_t, 4>, CosinePotential<float >>> read_local_potential(const toml::value& local);
#endif

} // mjolnir
#endif // MJOLNIR_READ_LOCAL_POTENTIAL_HPP
