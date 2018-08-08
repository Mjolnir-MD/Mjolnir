#define BOOST_TEST_MODULE "test_flexible_local_dihedral_potential"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/mjolnir/traits.hpp>
#include <mjolnir/potential/local/FlexibleLocalDihedralPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(FlexibleLocalDihedral_double)
{
    using real_type = double;
    constexpr static std::size_t N = 1000;
    constexpr static real_type   h = 1e-6;

    const real_type k  = 1.0;
    const std::array<real_type, 7> term{{2.2056, 0.2183, -0.0795, 0.0451,
                                        -0.3169, 0.0165, -0.1375}};

    mjolnir::FlexibleLocalDihedralPotential<real_type> flpd(k, term);

    const real_type x_min = -3.14;
    const real_type x_max = 3.14;
    const real_type dx = (x_max - x_min) / N;

    real_type x = x_min;
    for(std::size_t i=0; i<N; ++i)
    {
        const real_type pot1 = flpd.potential(x + h);
        const real_type pot2 = flpd.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = flpd.derivative(x);

        if(std::abs(deri) > h)
            BOOST_CHECK_CLOSE_FRACTION(dpot, deri, h);
        else
            BOOST_CHECK_SMALL(deri, h);

        x += dx;
    }
}

BOOST_AUTO_TEST_CASE(FlexibleLocalDihedral_float)
{
    using real_type = float;
    constexpr static std::size_t N   = 100;
    constexpr static real_type   h   = 1e-3;
    constexpr static real_type   tol = 5e-3;

    const real_type k  = 1.0;
    const std::array<real_type, 7> term{{2.2056, 0.2183, -0.0795, 0.0451,
                                        -0.3169, 0.0165, -0.1375}};

    mjolnir::FlexibleLocalDihedralPotential<real_type> flpd(k, term);

    const real_type x_min = -3.14;
    const real_type x_max = 3.14;
    const real_type dx = (x_max - x_min) / N;

    real_type x = x_min;
    for(std::size_t i=0; i<N; ++i)
    {
        const real_type pot1 = flpd.potential(x + h);
        const real_type pot2 = flpd.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = flpd.derivative(x);

        if(std::abs(deri) > h)
            BOOST_CHECK_CLOSE_FRACTION(dpot, deri, tol);
        else
            BOOST_CHECK_SMALL(deri, tol);
        x += dx;
    }
}
