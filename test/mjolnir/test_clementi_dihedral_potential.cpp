#define BOOST_TEST_MODULE "test_clementi_dihedral_potential"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/mjolnir/traits.hpp>
#include <mjolnir/potential/local/ClementiDihedralPotential.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(ClementiDihedral_double)
{
    using real_type = double;
    constexpr static std::size_t N  = 1000;
    constexpr static real_type   h  = 1e-6;
    constexpr static real_type   pi = mjolnir::math::constants<real_type>::pi;

    const real_type k1 = 1.0;
    const real_type k3 = 2.0;
    const real_type r0 = pi * 2. / 3.;

    mjolnir::ClementiDihedralPotential<real_type> c(k1, k3, r0);

    const real_type x_min = 0.;
    const real_type x_max = 2. * pi;
    const real_type dx = (x_max - x_min) / N;

    real_type x = x_min;
    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x_up = (x+h < 2. * pi) ? x+h : x+h - 2. * pi;
        const real_type x_bt = (x-h > 0) ? x-h : x-h + 2. * pi;
        const real_type pot1 = c.potential(x_up);
        const real_type pot2 = c.potential(x_bt);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = c.derivative(x);

        if(std::abs(deri) > h)
            BOOST_CHECK_CLOSE_FRACTION(dpot, deri, h);
        else
            BOOST_CHECK_SMALL(deri, h);

        x += dx;
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

    real_type x = x_min;
    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x_up = (x+h < 2.0f * pi) ? x+h : x+h - 2.0f * pi;
        const real_type x_bt = (x-h > 0.0f) ? x-h : x-h + 2.0f * pi;
        const real_type pot1 = c.potential(x_up);
        const real_type pot2 = c.potential(x_bt);
        const real_type dpot = (pot1 - pot2) / (2.0f * h);
        const real_type deri = c.derivative(x);

        if(std::abs(deri) > h)
        {
            BOOST_CHECK_CLOSE_FRACTION(dpot, deri, h);
        }
        else
        {
            BOOST_CHECK_SMALL(deri, h);
        }

        x += dx;
    }
}

