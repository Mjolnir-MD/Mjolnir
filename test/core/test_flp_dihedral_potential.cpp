#define BOOST_TEST_MODULE "test_flexible_local_dihedral_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/check_potential.hpp>
#include <mjolnir/forcefield/FLP/FlexibleLocalDihedralPotential.hpp>
#include <mjolnir/math/constants.hpp>

BOOST_AUTO_TEST_CASE(FlexibleLocalDihedral_double)
{
    using real_type = double;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-5;
    constexpr real_type tol = 1e-5;
    constexpr real_type   pi = mjolnir::math::constants<real_type>::pi();

    const real_type k  = 1.0;
    const std::array<real_type, 7> term{{
        2.2056, 0.2183, -0.0795, 0.0451, -0.3169, 0.0165, -0.1375
    }};

    mjolnir::FlexibleLocalDihedralPotential<real_type> pot(k, term);

    const real_type x_min = -pi;
    const real_type x_max =  pi;

    mjolnir::test::check_potential(pot, x_min, x_max, tol, h, N);
}

BOOST_AUTO_TEST_CASE(FlexibleLocalDihedral_float)
{
    using real_type = float;
    constexpr std::size_t N   = 100;
    constexpr real_type   h   = 1e-2;
    constexpr real_type   tol = 1e-2;
    constexpr real_type   pi = mjolnir::math::constants<real_type>::pi();

    const real_type k  = 1.0;
    const std::array<real_type, 7> term{{
        2.2056, 0.2183, -0.0795, 0.0451, -0.3169, 0.0165, -0.1375
    }};

    mjolnir::FlexibleLocalDihedralPotential<real_type> pot(k, term);

    const real_type x_min = -pi;
    const real_type x_max =  pi;

    mjolnir::test::check_potential(pot, x_min, x_max, tol, h, N);
}
