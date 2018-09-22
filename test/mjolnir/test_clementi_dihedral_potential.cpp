#define BOOST_TEST_MODULE "test_clementi_dihedral_potential"
#include <boost/test/included/unit_test.hpp>

#include <mjolnir/potential/local/ClementiDihedralPotential.hpp>
#include <mjolnir/math/constants.hpp>

BOOST_AUTO_TEST_CASE(ClementiDihedral_double)
{
    using real_type = double;
    constexpr std::size_t N  = 1000;
    constexpr real_type   h  = 1e-6;
    constexpr real_type   pi = mjolnir::math::constants<real_type>::pi;

    const real_type k1 = 1.0;
    const real_type k3 = 2.0;
    const real_type r0 = pi * 2. / 3.;

    mjolnir::ClementiDihedralPotential<real_type> c(k1, k3, r0);

    const real_type x_min = -pi;
    const real_type x_max =  pi;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + dx * i;
        const real_type x_up = (x+h < 2. * pi) ? x+h : x+h - 2 * pi;
        const real_type x_bt = (x-h > 0)       ? x-h : x-h + 2 * pi;
        const real_type pot1 = c.potential(x_up);
        const real_type pot2 = c.potential(x_bt);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = c.derivative(x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
        // periodicity
        BOOST_TEST(c.potential(x)  == c.potential (x + 2 * pi),
                   boost::test_tools::tolerance(1e-8));
        BOOST_TEST(c.derivative(x) == c.derivative(x + 2 * pi),
                   boost::test_tools::tolerance(1e-8));
    }
}

BOOST_AUTO_TEST_CASE(ClementiDihedral_float)
{
    using real_type = float;
    constexpr static std::size_t N  = 100;
    constexpr static real_type   h  = 1e-2f;
    constexpr static real_type   pi = mjolnir::math::constants<real_type>::pi;

    const real_type k1 = 1.0f;
    const real_type k3 = 2.0f;
    const real_type r0 = pi * 2. / 3.;

    mjolnir::ClementiDihedralPotential<real_type> c(k1, k3, r0);

    const real_type x_min = 0.f;
    const real_type x_max = 2.f * pi;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + dx * i;
        const real_type x_up = (x+h < 2.0f * pi) ? x+h : x+h - 2 * pi;
        const real_type x_bt = (x-h > 0.0f)      ? x-h : x-h + 2 * pi;
        const real_type pot1 = c.potential(x_up);
        const real_type pot2 = c.potential(x_bt);
        const real_type dpot = (pot1 - pot2) / (2.0f * h);
        const real_type deri = c.derivative(x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
        // periodicity
        BOOST_TEST(c.potential(x)  == c.potential (x + 2 * pi),
                   boost::test_tools::tolerance(h));
        BOOST_TEST(c.derivative(x) == c.derivative(x + 2 * pi),
                   boost::test_tools::tolerance(h));
   }
}

