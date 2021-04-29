#define BOOST_TEST_MODULE "test_lennard_jones_attractive_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/forcefield/global/LennardJonesPotential.hpp>
#include <mjolnir/forcefield/global/LennardJonesAttractivePotential.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

BOOST_AUTO_TEST_CASE(LennardJones_double)
{
    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = typename traits_type::real_type;
    using molecule_id_type = mjolnir::Topology::molecule_id_type;
    using group_id_type    = mjolnir::Topology::group_id_type;
    using potential_type   = mjolnir::LennardJonesAttractivePotential<traits_type>;
    using parameter_type   = potential_type::parameter_type;

    constexpr static std::size_t N = 10000;
    constexpr static real_type   h = 1e-6;

    const real_type sigma   = 3.0;
    const real_type epsilon = 1.0;
    const parameter_type param{sigma, epsilon};
    potential_type lj{
        potential_type::default_cutoff(),
        {{0, param}, {1, param}}, {},
        mjolnir::IgnoreMolecule<molecule_id_type>("Nothing"),
        mjolnir::IgnoreGroup   <group_id_type   >({})
    };

    const real_type x_min = std::pow(2.0, 1.0/6.0) * sigma;
    const real_type x_max = lj.cutoff_ratio() * sigma;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=1; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = lj.potential(0, 1, x + h);
        const real_type pot2 = lj.potential(0, 1, x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = lj.derivative(0, 1, x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
    }

    // -epsilon if x < 2^1/6 sigma
    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x = (double(i) / N) * x_min;
        const real_type E = lj.potential(0, 1, x);
        BOOST_TEST(E == -epsilon, boost::test_tools::tolerance(h));
    }

    mjolnir::LennardJonesPotential<traits_type> ref{
        potential_type::default_cutoff(),
        {{0, param}, {1, param}}, {},
        mjolnir::IgnoreMolecule<molecule_id_type>("Nothing"),
        mjolnir::IgnoreGroup   <group_id_type   >({})
    };

    for(std::size_t i=1; i<N; ++i)
    {
        const real_type x      = x_min + i * dx;
        const real_type pot    = lj.potential(0, 1, x);
        const real_type ref_lj = ref.potential(0, 1, x);
        BOOST_TEST(pot == ref_lj, boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(LennardJones_float)
{
    using traits_type = mjolnir::SimulatorTraits<float, mjolnir::UnlimitedBoundary>;
    using real_type   = typename traits_type::real_type;
    using molecule_id_type = mjolnir::Topology::molecule_id_type;
    using group_id_type    = mjolnir::Topology::group_id_type;
    using potential_type   = mjolnir::LennardJonesAttractivePotential<traits_type>;
    using parameter_type   = potential_type::parameter_type;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 0.002f;
    constexpr real_type tol = 0.005f;

    const real_type sigma   = 3.0f;
    const real_type epsilon = 1.0f;
    const parameter_type param{sigma, epsilon};
    potential_type lj{potential_type::default_cutoff(),
        {{0, param}, {1, param}}, {},
        mjolnir::IgnoreMolecule<molecule_id_type>("Nothing"),
        mjolnir::IgnoreGroup   <group_id_type   >({})
    };
    const real_type cutoff = lj.cutoff_ratio();

    const real_type x_min = std::pow(2.0f, 1.0f/6.0f) * sigma;
    const real_type x_max = cutoff * sigma;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=1; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = lj.potential(0, 1, x + h);
        const real_type pot2 = lj.potential(0, 1, x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = lj.derivative(0, 1, x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(tol));
    }
}


