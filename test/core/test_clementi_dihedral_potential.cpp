#define BOOST_TEST_MODULE "test_clementi_dihedral_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/check_potential.hpp>
#include <mjolnir/forcefield/local/ClementiDihedralPotential.hpp>
#include <mjolnir/math/constants.hpp>

BOOST_AUTO_TEST_CASE(ClementiDihedral_double)
{
    using real_type = double;
    constexpr std::size_t N  = 1000;
    constexpr real_type   h  = 1e-6;
    constexpr real_type tol  = 1e-6;
    constexpr real_type   pi = mjolnir::math::constants<real_type>::pi();
    constexpr real_type   two_pi = 2.0 * pi;

    const real_type k1 = 1.0;
    const real_type k3 = 2.0;
    const real_type r0 = pi * 2. / 3.;

    mjolnir::ClementiDihedralPotential<real_type> potential(k1, k3, r0);

    const real_type x_min = -pi;
    const real_type x_max =  pi;

    mjolnir::test::check_potential(potential, x_min, x_max, tol, h, N);

    // periodicity
    const real_type dx = (x_max - x_min) / N;
    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x = x_min + dx * i;
        BOOST_TEST(potential.potential(x) == potential.potential(x + two_pi),
                   boost::test_tools::tolerance(1e-8));
    }
}

BOOST_AUTO_TEST_CASE(ClementiDihedral_float)
{
    using real_type = float;
    constexpr std::size_t N  = 100;
    constexpr real_type   h  = 1e-2f;
    constexpr real_type tol  = 1e-2f;
    constexpr real_type   pi = mjolnir::math::constants<real_type>::pi();
    constexpr real_type   two_pi = 2.0 * pi;

    const real_type k1 = 1.0f;
    const real_type k3 = 2.0f;
    const real_type r0 = pi * 2. / 3.;

    mjolnir::ClementiDihedralPotential<real_type> potential(k1, k3, r0);

    const real_type x_min = 0.f;
    const real_type x_max = 2.f * pi;

    mjolnir::test::check_potential(potential, x_min, x_max, tol, h, N);

    // periodicity
    const real_type dx = (x_max - x_min) / N;
    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + dx * i;
        BOOST_TEST(potential.potential(x) == potential.potential(x + two_pi),
                   boost::test_tools::tolerance(h));
   }
}

