#define BOOST_TEST_MODULE "test_go1012_potential"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/potential/local/Go1012ContactPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(GoContact_double)
{
    using real_type = double;
    constexpr static std::size_t N = 1000;
    constexpr static real_type   h = 1e-6;

    const real_type e  = 1.0;
    const real_type r0 = 5.0;

    mjolnir::Go1012ContactPotential<real_type> Go(e, r0);

    const real_type x_min = 0.8 * r0;
    const real_type x_max = 5.0 * r0;
    const real_type dx = (x_max - x_min) / N;

    real_type x = x_min;
    for(std::size_t i=0; i<N; ++i)
    {
        const real_type pot1 = Go.potential(x + h);
        const real_type pot2 = Go.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = Go.derivative(x);

        if(std::abs(deri) > h)
            BOOST_CHECK_CLOSE_FRACTION(dpot, deri, h);
        else
            BOOST_CHECK_SMALL(deri, h);
        x += dx;
    }
}

BOOST_AUTO_TEST_CASE(GoContact_float)
{
    using real_type = double;
    constexpr static std::size_t       N    = 100;
    constexpr static real_type h    = 1e-3;

    const real_type e  = 1.0;
    const real_type r0 = 5.0;

    mjolnir::Go1012ContactPotential<real_type> Go(e, r0);

    const real_type x_min = 0.8 * r0;
    const real_type x_max = 5.0 * r0;
    const real_type dx = (x_max - x_min) / N;

    real_type x = x_min;
    for(std::size_t i=0; i<N; ++i)
    {
        const real_type pot1 = Go.potential(x + h);
        const real_type pot2 = Go.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = Go.derivative(x);

        if(std::abs(deri) > h)
            BOOST_CHECK_CLOSE_FRACTION(dpot, deri, h);
        else
            BOOST_CHECK_SMALL(deri, h);
        x += dx;
    }
}
