#define BOOST_TEST_MODULE "test_uniform_cubic_pan_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/check_potential.hpp>
#include <mjolnir/forcefield/global/UniformCubicPanPotential.hpp>

BOOST_AUTO_TEST_CASE(UniformCubicPan_double)
{
    using real_type   = double;
    using potential_type = mjolnir::UniformCubicPanPotential<real_type>;
    using parameter_type = potential_type::parameter_type;

    constexpr std::size_t N = 10000;
    constexpr real_type   h = 1e-6;
    constexpr real_type tol = 1e-6;

    constexpr real_type epsilon = 1.0;
    constexpr real_type v0      = 5.0;
    constexpr real_type range   = 5.0;

    potential_type potential(epsilon, v0, range);

    const real_type x_min = 0.8 * v0;
    const real_type x_max = (v0 + range) * 1.2;

    mjolnir::test::check_potential(potential, parameter_type{},
                                   x_min, x_max, tol, h, N);
}

BOOST_AUTO_TEST_CASE(UniformCubicPan_float)
{
    using real_type   = float;
    using potential_type = mjolnir::UniformCubicPanPotential<real_type>;
    using parameter_type = potential_type::parameter_type;

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 0.002f;
    constexpr real_type tol = 0.005f;

    constexpr real_type epsilon = 1.0;
    constexpr real_type v0      = 5.0;
    constexpr real_type range   = 5.0;

    potential_type potential(epsilon, v0, range);

    const real_type x_min = 0.8 * v0;
    const real_type x_max = (v0 + range) * 1.2;

    mjolnir::test::check_potential(potential, parameter_type{},
                                   x_min, x_max, tol, h, N);
}
