#define BOOST_TEST_MODULE "test_uniform_lennard_jones_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/forcefield/global/UniformLennardJonesPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(UniformLennardJones_double)
{
    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = typename traits_type::real_type;
    using molecule_id_type = mjolnir::Topology::molecule_id_type;
    using group_id_type    = mjolnir::Topology::group_id_type;
    constexpr std::size_t N = 10000;
    constexpr real_type   h = 1e-6;

    constexpr real_type sigma   = 3.0;
    constexpr real_type epsilon = 1.0;
    mjolnir::UniformLennardJonesPotential<traits_type> lj{
        sigma, epsilon,
        mjolnir::UniformLennardJonesPotential<traits_type>::default_cutoff(), {}, {},
        mjolnir::IgnoreMolecule<molecule_id_type>("Nothing"),
        mjolnir::IgnoreGroup   <group_id_type   >({})
    };
    const real_type cutoff = lj.cutoff_ratio();

    const real_type x_min = 0.8 * sigma;
    const real_type x_max = cutoff * sigma;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = lj.potential(0, 1, x + h);
        const real_type pot2 = lj.potential(0, 1, x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = lj.derivative(0, 1, x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(UniformLennardJones_float)
{
    using traits_type = mjolnir::SimulatorTraits<float, mjolnir::UnlimitedBoundary>;
    using real_type   = typename traits_type::real_type;
    using molecule_id_type = mjolnir::Topology::molecule_id_type;
    using group_id_type    = mjolnir::Topology::group_id_type;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 0.002;
    constexpr real_type tol = 0.005;

    constexpr real_type sigma   = 3.0;
    constexpr real_type epsilon = 1.0;
    mjolnir::UniformLennardJonesPotential<traits_type> lj{
        sigma, epsilon,
        mjolnir::UniformLennardJonesPotential<traits_type>::default_cutoff(), {}, {},
        mjolnir::IgnoreMolecule<molecule_id_type>("Nothing"),
        mjolnir::IgnoreGroup   <group_id_type   >({})
    };
    const real_type cutoff = lj.cutoff_ratio();

    const real_type x_min = 0.8 * sigma;
    const real_type x_max = cutoff * sigma;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = lj.potential(0, 1, x + h);
        const real_type pot2 = lj.potential(0, 1, x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = lj.derivative(0, 1, x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(tol));
    }
}
