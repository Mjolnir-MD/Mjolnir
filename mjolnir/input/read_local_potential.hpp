#ifndef MJOLNIR_READ_LOCAL_POTENTIAL_HPP
#define MJOLNIR_READ_LOCAL_POTENTIAL_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/potential/local/HarmonicPotential.hpp>
#include <mjolnir/potential/local/Go1012ContactPotential.hpp>
#include <mjolnir/potential/local/ClementiDihedralPotential.hpp>
#include <mjolnir/potential/local/GaussianPotential.hpp>
#include <mjolnir/potential/local/PeriodicGaussianPotential.hpp>
#include <mjolnir/potential/local/FlexibleLocalAnglePotential.hpp>
#include <mjolnir/potential/local/FlexibleLocalDihedralPotential.hpp>
#include <mjolnir/potential/local/SumLocalPotential.hpp>
#include <mjolnir/core/Topology.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// this reads...
// parameters = [
// {indices = [1, 2], k = 10.0, v0 = 1.0} <- a portion of this table, k and v0.
// ]
template<typename realT>
HarmonicPotential<realT> read_harmonic_potential(const toml::value& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    const auto k  = toml::find<real_type>(param, "k" );
    const auto v0 = toml::find<real_type>(param, "v0");

    MJOLNIR_LOG_INFO("HarmonicPotential = {v0 = ", v0, ", k = ", k, '}');
    return HarmonicPotential<realT>(k, v0);
}

template<typename realT>
Go1012ContactPotential<realT>
read_go1012_contact_potential(const toml::value& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    const auto k  = toml::find<real_type>(param, "k");
    const auto v0 = toml::find<real_type>(param, "v0");

    MJOLNIR_LOG_INFO("GoContactPotential = {v0 = ", v0, ", k = ", k, '}');
    return Go1012ContactPotential<realT>(k, v0);
}

template<typename realT>
GaussianPotential<realT> read_gaussian_potential(const toml::value& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    const auto v0    = toml::find<real_type>(param, "v0");
    const auto k     = toml::find<real_type>(param, "k");
    const auto sigma = toml::expect<real_type>(param, u8"σ").or_other(
                       toml::expect<real_type>(param, "sigma")).unwrap();
    MJOLNIR_LOG_INFO("GaussianPotential = {v0 = ", v0, ", k = ", k,
                     ", sigma = ", sigma, '}');
    return GaussianPotential<realT>(k, sigma, v0);
}

template<typename realT>
PeriodicGaussianPotential<realT>
read_periodic_gaussian_potential(const toml::value& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    const auto v0    = toml::find<real_type>(param, "v0");
    const auto k     = toml::find<real_type>(param, "k");
    const auto sigma = toml::expect<real_type>(param, u8"σ").or_other(
                       toml::expect<real_type>(param, "sigma")).unwrap();

    MJOLNIR_LOG_INFO("PeriodicGaussianPotential = {v0 = ", v0, ", k = ", k,
                     ", sigma = ", sigma, '}');
    return PeriodicGaussianPotential<realT>(k, sigma, v0);
}

template<typename realT>
FlexibleLocalAnglePotential<realT>
read_flexible_local_angle_potential(const toml::value& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    const auto k     = toml::find<real_type                >(param, "k");
    const auto term1 = toml::find<std::array<real_type, 10>>(param, "y");
    const auto term2 = toml::find<std::array<real_type, 10>>(param, "d2y");

    MJOLNIR_LOG_INFO("FlexibleLocalAngle = {k = ", k,
                     ", y = ", term1, ", d2y = ", term2, '}');
    return FlexibleLocalAnglePotential<realT>(k, term1, term2);
}

template<typename realT>
ClementiDihedralPotential<realT>
read_clementi_dihedral_potential(const toml::value& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    const auto v0 = toml::find<real_type>(param, "v0");
    const auto k1 = toml::find<real_type>(param, "k1");
    const auto k3 = toml::find<real_type>(param, "k3");

    MJOLNIR_LOG_INFO("ClementiDihedral = {v0 = ", v0,
                     ", k1 = ", k1, ", k3 = ", k3, '}');
    return ClementiDihedralPotential<realT>(k1, k3, v0);
}

template<typename realT>
FlexibleLocalDihedralPotential<realT>
read_flexible_local_dihedral_potential(const toml::value& param)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = realT;
    auto k    = toml::find<real_type               >(param, "k");
    auto term = toml::find<std::array<real_type, 7>>(param, "coef");

    MJOLNIR_LOG_INFO("FlexibleLocalDihedral = {k = ", k, ", coef = ", term, '}');
    return FlexibleLocalDihedralPotential<realT>(k, term);
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
    static HarmonicPotential<realT> invoke(const toml::value& param)
    {
        return read_harmonic_potential<realT>(param);
    }
};
template<typename realT>
struct read_local_potential_impl<Go1012ContactPotential<realT>>
{
    static Go1012ContactPotential<realT> invoke(const toml::value& param)
    {
        return read_go1012_contact_potential<realT>(param);
    }
};
template<typename realT>
struct read_local_potential_impl<GaussianPotential<realT>>
{
    static GaussianPotential<realT> invoke(const toml::value& param)
    {
        return read_gaussian_potential<realT>(param);
    }
};
template<typename realT>
struct read_local_potential_impl<PeriodicGaussianPotential<realT>>
{
    static PeriodicGaussianPotential<realT> invoke(const toml::value& param)
    {
        return read_periodic_gaussian_potential<realT>(param);
    }
};
template<typename realT>
struct read_local_potential_impl<FlexibleLocalAnglePotential<realT>>
{
    static FlexibleLocalAnglePotential<realT> invoke(const toml::value& param)
    {
        return read_flexible_local_angle_potential<realT>(param);
    }
};
template<typename realT>
struct read_local_potential_impl<ClementiDihedralPotential<realT>>
{
    static ClementiDihedralPotential<realT> invoke(const toml::value& param)
    {
        return read_clementi_dihedral_potential<realT>(param);
    }
};
template<typename realT>
struct read_local_potential_impl<FlexibleLocalDihedralPotential<realT>>
{
    static FlexibleLocalDihedralPotential<realT> invoke(const toml::value& param)
    {
        return read_flexible_local_dihedral_potential<realT>(param);
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
    MJOLNIR_SCOPE(read_local_potential(), 0);
    MJOLNIR_LOG_INFO("as ", N, "-body interaction");

    using indices_t                = std::array<std::size_t, N>;
    using indices_potential_pair_t = std::pair<indices_t, potentialT>;

    const auto& params = toml::find<toml::array>(local, "parameters");
    MJOLNIR_LOG_NOTICE("-- ", params.size(), " interactions are found.");

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());
    for(const auto& item : params)
    {
        const auto indices = toml::find<indices_t>(item, "indices");
        MJOLNIR_LOG_INFO_NO_LF("idxs = ", indices, ", ");

        retval.emplace_back(indices,
            detail::read_local_potential_impl<potentialT>::invoke(item));
    }
    return retval;
}

template<std::size_t N, typename realT,
         typename potentialT1, typename potentialT2>
std::vector<std::pair<std::array<std::size_t, N>,
            SumLocalPotential<realT, potentialT1, potentialT2>>>
read_local_potentials(const toml::value& local,
                      const std::string& p1, const std::string& p2)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE(read_local_potentials(), 0);
    MJOLNIR_LOG_INFO("as ", N, "-body interaction");

    using indices_t                = std::array<std::size_t, N>;
    using indices_potential_pair_t = std::pair<indices_t,
        SumLocalPotential<realT, potentialT1, potentialT2>>;

    const auto& params = toml::find<toml::Array>(local, "parameters");
    MJOLNIR_LOG_NOTICE("-- ", params.size(), " interactions are found.");

    std::vector<indices_potential_pair_t> retval;
    retval.reserve(params.size());
    for(const auto& item : params)
    {
        const auto indices = toml::find<indices_t>(item, "indices");
        MJOLNIR_LOG_INFO("idxs = ", indices);

        const auto& pot1 = toml::find<toml::value>(item, p1);
        const auto& pot2 = toml::find<toml::value>(item, p2);

        retval.emplace_back(indices,
            SumLocalPotential<realT, potentialT1, potentialT2>(
                detail::read_local_potential_impl<potentialT1>::invoke(pot1),
                detail::read_local_potential_impl<potentialT2>::invoke(pot2)
            ));
    }
    return retval;
}

} // mjolnir
#endif // MJOLNIR_READ_LOCAL_POTENTIAL_HPP
