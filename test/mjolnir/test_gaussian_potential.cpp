#define BOOST_TEST_MODULE "test_gaussian_potential"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/mjolnir/traits.hpp>
#include <mjolnir/potential/local/GaussianPotential.hpp>
#include <mjolnir/util/make_unique.hpp>


BOOST_AUTO_TEST_CASE(Gaussian_double)
{
    typedef mjolnir::test::traits<double> traits;
    constexpr static std::size_t       N    = 1000;
    constexpr static traits::real_type h    = 1e-6;
    const traits::real_type e  = 2.0;
    const traits::real_type w  = 0.15;
    const traits::real_type r0 = 7.0;

    mjolnir::GaussianPotential<traits> gaussian(e, w, r0);

    const traits::real_type x_min = 0.5 * r0;
    const traits::real_type x_max = 1.5 * r0;
    const traits::real_type dx = (x_max - x_min) / N;

    traits::real_type x = x_min;
    for(std::size_t i=0; i<N; ++i)
    {
        const traits::real_type pot1 = gaussian.potential(x + h);
        const traits::real_type pot2 = gaussian.potential(x - h);
        const traits::real_type dpot = (pot1 - pot2) / (2 * h);
        const traits::real_type deri = gaussian.derivative(x);

        if(std::abs(deri) > h)
            BOOST_CHECK_CLOSE_FRACTION(dpot, deri, h);
        else
            BOOST_CHECK_SMALL(deri, h);
        x += dx;
    }
}

BOOST_AUTO_TEST_CASE(Gaussian_float)
{
    typedef mjolnir::test::traits<float> traits;
    constexpr static std::size_t       N    = 100;
    constexpr static traits::real_type h    = 1e-3;
    const traits::real_type e  = 2.0;
    const traits::real_type w  = 0.15;
    const traits::real_type r0 = 7.0;

    mjolnir::GaussianPotential<traits> gaussian(e, w, r0);

    const traits::real_type x_min = 0.5 * r0;
    const traits::real_type x_max = 1.5 * r0;
    const traits::real_type dx = (x_max - x_min) / N;

    traits::real_type x = x_min;
    for(std::size_t i=0; i<N; ++i)
    {
        const traits::real_type pot1 = gaussian.potential(x + h);
        const traits::real_type pot2 = gaussian.potential(x - h);
        const traits::real_type dpot = (pot1 - pot2) / (2 * h);
        const traits::real_type deri = gaussian.derivative(x);

        if(std::abs(deri) > h)
            BOOST_CHECK_CLOSE_FRACTION(dpot, deri, h);
        else
            BOOST_CHECK_SMALL(deri, h);
        x += dx;
    }
}
