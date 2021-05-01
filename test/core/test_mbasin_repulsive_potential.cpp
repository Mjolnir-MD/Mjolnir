#define BOOST_TEST_MODULE "test_mbasin_repulsive_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/forcefield/local/GoContactPotential.hpp>
#include <mjolnir/forcefield/MultipleBasin/MBasinRepulsivePotential.hpp>
#include <iomanip>

BOOST_AUTO_TEST_CASE(MBasinrepulsive)
{
    using real_type = double;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;

    const real_type e  = 1.0;
    const real_type r0 = 5.0;

    mjolnir::GoContactPotential<real_type>       ref(e, r0);
    mjolnir::MBasinRepulsivePotential<real_type> pot(e, r0);

    {
        const real_type x_min = r0 * 4.0 / 6.0;
        const real_type x_max = r0 * std::sqrt(5.0 / 6.0);
        const real_type dx = (x_max - x_min) / N;
        for(std::size_t i=0; i<N; ++i)
        {
            const real_type x = x_min + i * dx;
            const real_type p1 = pot.potential(x);
            const real_type p2 = ref.potential(x);
            const real_type d1 = pot.derivative(x);
            const real_type d2 = ref.derivative(x);

            BOOST_TEST(p1 == p2, boost::test_tools::tolerance(h));
            BOOST_TEST(d1 == d2, boost::test_tools::tolerance(h));
        }
    }
    {
        const real_type x_min = r0 * std::sqrt(5.0 / 6.0);
        const real_type x_max = r0 * 2.5;
        const real_type dx = (x_max - x_min) / N;
        for(std::size_t i=1; i<N; ++i)
        {
            const real_type x = x_min + i * dx;
            const real_type p = pot.potential(x);
            const real_type d = pot.derivative(x);

            BOOST_TEST(p == 0.0, boost::test_tools::tolerance(h));
            BOOST_TEST(d == 0.0, boost::test_tools::tolerance(h));
        }
    }
}
