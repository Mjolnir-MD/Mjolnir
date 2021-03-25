#define BOOST_TEST_MODULE "test_inverse_power_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/forcefield/global/InversePowerPotential.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

BOOST_AUTO_TEST_CASE(InversePower_double)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_inverse_power_potential.log");

    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = typename traits_type::real_type;
    using integer_type = typename mjolnir::InversePowerPotential<traits_type>::integer_type;
    using molecule_id_type = mjolnir::Topology::molecule_id_type;
    using group_id_type    = mjolnir::Topology::group_id_type;
    constexpr std::size_t N = 10000;
    constexpr real_type   h = 1e-6;

    const real_type sigma   = 3.0;
    const real_type epsilon = 1.0;
    const integer_type n    = 5;
    mjolnir::InversePowerPotential<traits_type> inp{
        epsilon, n, mjolnir::InversePowerPotential<traits_type>::default_cutoff(n),
        {{0, sigma}, {1, sigma}}, {},
        mjolnir::IgnoreMolecule<molecule_id_type>("Nothing"),
        mjolnir::IgnoreGroup   <group_id_type   >({})
    };

    const real_type x_min = 0.8 * sigma;
    const real_type x_max = inp.cutoff_ratio() * sigma;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = inp.potential(0, 1, x + h);
        const real_type pot2 = inp.potential(0, 1, x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = inp.derivative(0, 1, x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(InversePower_float)
{
    using traits_type = mjolnir::SimulatorTraits<float, mjolnir::UnlimitedBoundary>;
    using real_type   = typename traits_type::real_type;
    using integer_type = typename mjolnir::InversePowerPotential<traits_type>::integer_type;
    using molecule_id_type = mjolnir::Topology::molecule_id_type;
    using group_id_type    = mjolnir::Topology::group_id_type;
    constexpr static std::size_t N = 1000;
    constexpr static real_type   h = 0.002;
    constexpr static real_type tol = 0.005;

    const real_type sigma   = 3.0;
    const real_type epsilon = 1.0;
    const integer_type n    = 5;
    mjolnir::InversePowerPotential<traits_type> inp{
        epsilon, n, mjolnir::InversePowerPotential<traits_type>::default_cutoff(n),
        {{0, sigma}, {1, sigma}}, {},
        mjolnir::IgnoreMolecule<molecule_id_type>("Nothing"),
        mjolnir::IgnoreGroup   <group_id_type   >({})
    };
    const real_type cutoff = inp.cutoff_ratio();

    const real_type x_min = 0.8f   * sigma;
    const real_type x_max = cutoff * sigma;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = inp.potential(0, 1, x + h);
        const real_type pot2 = inp.potential(0, 1, x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = inp.derivative(0, 1, x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(tol));
    }
}

