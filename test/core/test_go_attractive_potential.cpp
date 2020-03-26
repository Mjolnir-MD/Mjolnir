#define BOOST_TEST_MODULE "test_attractive_gocontact_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/forcefield/local/GoContactAttractivePotential.hpp>
#include <iomanip>

BOOST_AUTO_TEST_CASE(GoContactAttractive_double)
{
    using real_type = double;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;

    const real_type e  = 1.0;
    const real_type r0 = 5.0;

    mjolnir::GoContactAttractivePotential<real_type> pot(e, r0);

    {
        const real_type x_min = 0.0;
        const real_type x_max = r0;
        const real_type dx    = (x_max - x_min) / N;

        for(std::size_t i=0; i<N; ++i)
        {
            const real_type x = x_min + i * dx;
            const real_type p = pot.potential(x);
            const real_type d = pot.derivative(x);
            BOOST_TEST(p == -e);
            BOOST_TEST(d == 0.0);
        }
    }
    {
        const real_type x_min = r0 * 1.01;
        const real_type x_max = pot.cutoff();
        const real_type dx = (x_max - x_min) / N;
        for(std::size_t i=0; i<N; ++i)
        {
            const real_type x    = x_min + i * dx;
            const real_type pot1 = pot.potential(x + h);
            const real_type pot2 = pot.potential(x - h);
            const real_type dpot = (pot1 - pot2) / (2 * h);
            const real_type deri = pot.derivative(x);

            BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
        }
    }
    {
        const real_type x_min = pot.cutoff() + 1.01;
        const real_type x_max = 2 * pot.cutoff();
        const real_type dx = (x_max - x_min) / N;
        for(std::size_t i=0; i<N; ++i)
        {
            const real_type x = x_min + i * dx;
            BOOST_TEST(pot.potential(x)  == 0.0, boost::test_tools::tolerance(h));
            BOOST_TEST(pot.derivative(x) == 0.0, boost::test_tools::tolerance(h));
        }
    }
}

BOOST_AUTO_TEST_CASE(GoContactAttractive_float)
{
    using real_type = float;
    constexpr std::size_t N = 100;
    constexpr real_type   h = 1e-3;

    const real_type e  = 1.0;
    const real_type r0 = 5.0;

    mjolnir::GoContactAttractivePotential<real_type> pot(e, r0);

    {
        const real_type x_min = 0.0;
        const real_type x_max = r0;
        const real_type dx    = (x_max - x_min) / N;

        for(std::size_t i=0; i<N; ++i)
        {
            const real_type x = x_min + i * dx;
            const real_type p = pot.potential(x);
            const real_type d = pot.derivative(x);
            BOOST_TEST(p == -e);
            BOOST_TEST(d == 0.0);
        }
    }
    {
        const real_type x_min = r0 * 1.01f;
        const real_type x_max = pot.cutoff();
        const real_type dx = (x_max - x_min) / N;
        for(std::size_t i=0; i<N; ++i)
        {
            const real_type x    = x_min + i * dx;
            const real_type pot1 = pot.potential(x + h);
            const real_type pot2 = pot.potential(x - h);
            const real_type dpot = (pot1 - pot2) / (2 * h);
            const real_type deri = pot.derivative(x);

            BOOST_TEST(dpot == deri, boost::test_tools::tolerance(1e-2f));
        }
    }
    {
        const real_type x_min = pot.cutoff() + 1.01f;
        const real_type x_max = 2 * pot.cutoff();
        const real_type dx = (x_max - x_min) / N;
        for(std::size_t i=0; i<N; ++i)
        {
            const real_type x = x_min + i * dx;
            BOOST_TEST(pot.potential(x)  == 0.0, boost::test_tools::tolerance(h));
            BOOST_TEST(pot.derivative(x) == 0.0, boost::test_tools::tolerance(h));
        }
    }
}
