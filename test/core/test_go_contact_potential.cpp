#define BOOST_TEST_MODULE "test_gocontact_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/forcefield/local/GoContactPotential.hpp>

BOOST_AUTO_TEST_CASE(GoContact_double)
{
    using real_type = double;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;

    const real_type e  = 1.0;
    const real_type r0 = 5.0;

    mjolnir::GoContactPotential<real_type> Go(e, r0);

    const real_type x_min = 0.8 * r0;
    const real_type x_max = 5.0 * r0;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = Go.potential(x + h);
        const real_type pot2 = Go.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = Go.derivative(x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(GoContact_float)
{
    using real_type = double;
    constexpr std::size_t N = 100;
    constexpr real_type   h = 1e-3;

    const real_type e  = 1.0;
    const real_type r0 = 5.0;

    mjolnir::GoContactPotential<real_type> Go(e, r0);

    const real_type x_min = 0.8 * r0;
    const real_type x_max = 5.0 * r0;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = Go.potential(x + h);
        const real_type pot2 = Go.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = Go.derivative(x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
    }
}
