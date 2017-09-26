#define BOOST_TEST_MODULE "test_flexible_local_angle_potential"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/potential/local/FlexibleLocalAnglePotential.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(FlexibleLocalAngle_double)
{
    typedef mjolnir::SimulatorTraitsBase<double, mjolnir::UnlimitedBoundary> traits;
    constexpr static std::size_t       N = 1000;
    constexpr static traits::real_type h = 1e-6;

    typedef typename traits::real_type real_type;
    const real_type k  = 1.0;
    const std::array<real_type, 10> term1{{5.00, 1.34, 0.84, 1.17, 0.82, 1.00,
                                           1.27, 1.52, 3.20, 10.00}};
    const std::array<real_type, 10> term2{{0.00, 151.96, 14.61, -46.89, 39.04,
                                          -4.86, -1.86, 8.38, 250.03, 0.00}};

    mjolnir::FlexibleLocalAnglePotential<traits> flpa(k, term1, term2);

    const real_type x_min = 1.309;
    const real_type x_max = 2.87979;
    const real_type dx = (x_max - x_min) / N;

    real_type x = x_min + dx;
    for(std::size_t i=1; i<N; ++i)
    {
        const real_type pot1 = flpa.potential(x + h);
        const real_type pot2 = flpa.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = flpa.derivative(x);

        if(std::abs(deri) > h)
            BOOST_CHECK_CLOSE_FRACTION(dpot, deri, h);
        else
            BOOST_CHECK_SMALL(deri, h);

        x += dx;
    }
}

BOOST_AUTO_TEST_CASE(FlexibleLocalAngle_float)
{
    typedef mjolnir::SimulatorTraitsBase<float, mjolnir::UnlimitedBoundary> traits;
    constexpr static std::size_t       N = 100;
    constexpr static traits::real_type h = 1e-3;

    typedef typename traits::real_type real_type;
    const real_type k  = 1.0;
    const std::array<real_type, 10> term1{{5.00, 1.34, 0.84, 1.17, 0.82, 1.00,
                                           1.27, 1.52, 3.20, 10.00}};
    const std::array<real_type, 10> term2{{0.00, 151.96, 14.61, -46.89, 39.04,
                                          -4.86, -1.86, 8.38, 250.03, 0.00}};

    mjolnir::FlexibleLocalAnglePotential<traits> flpa(k, term1, term2);

    const real_type x_min = 1.309;
    const real_type x_max = 2.87979;
    const real_type dx = (x_max - x_min) / N;

    real_type x = x_min + dx;
    for(std::size_t i=1; i<N; ++i)
    {
        const real_type pot1 = flpa.potential(x + h);
        const real_type pot2 = flpa.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = flpa.derivative(x);

        if(std::abs(deri) > h)
            BOOST_CHECK_CLOSE_FRACTION(dpot, deri, h);
        else
            BOOST_CHECK_SMALL(deri, h);

        x += dx;
    }
}
