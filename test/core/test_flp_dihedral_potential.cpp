#define BOOST_TEST_MODULE "test_flexible_local_dihedral_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/potential/local/FlexibleLocalDihedralPotential.hpp>
#include <mjolnir/math/constants.hpp>

BOOST_AUTO_TEST_CASE(FlexibleLocalDihedral_double)
{
    using real_type = double;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr real_type   pi = mjolnir::math::constants<real_type>::pi();

    const real_type k  = 1.0;
    const std::array<real_type, 7> term{{
        2.2056, 0.2183, -0.0795, 0.0451, -0.3169, 0.0165, -0.1375
    }};

    mjolnir::FlexibleLocalDihedralPotential<real_type> flpd(k, term);

    const real_type x_min = -pi;
    const real_type x_max =  pi;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + dx * i;
        const real_type pot1 = flpd.potential(x + h);
        const real_type pot2 = flpd.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = flpd.derivative(x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(1e-5));
        BOOST_TEST(flpd.potential (x) == flpd.potential (x + 2 * pi),
                   boost::test_tools::tolerance(h));
        BOOST_TEST(flpd.derivative(x) == flpd.derivative(x + 2 * pi),
                   boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(FlexibleLocalDihedral_float)
{
    using real_type = float;
    constexpr std::size_t N   = 100;
    constexpr real_type   h   = 1e-3;
    constexpr real_type   tol = 5e-2;
    constexpr real_type   pi = mjolnir::math::constants<real_type>::pi();

    const real_type k  = 1.0;
    const std::array<real_type, 7> term{{
        2.2056, 0.2183, -0.0795, 0.0451, -0.3169, 0.0165, -0.1375
    }};

    mjolnir::FlexibleLocalDihedralPotential<real_type> flpd(k, term);

    const real_type x_min = -pi;
    const real_type x_max =  pi;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + dx * i;
        const real_type pot1 = flpd.potential(x + h);
        const real_type pot2 = flpd.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = flpd.derivative(x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(tol));
        BOOST_TEST(flpd.potential (x) == flpd.potential (x + 2 * pi),
                   boost::test_tools::tolerance(tol));
        BOOST_TEST(flpd.derivative(x) == flpd.derivative(x + 2 * pi),
                   boost::test_tools::tolerance(tol));
    }
}
