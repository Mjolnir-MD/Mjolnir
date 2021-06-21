#define BOOST_TEST_MODULE "test_isolf_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>
#include <mjolnir/forcefield/iSoLF/iSoLFAttractivePotential.hpp>

BOOST_AUTO_TEST_CASE(iSoLF_double)
{
    using real_type = double;
    using potential_type = mjolnir::iSoLFAttractivePotential<real_type>;
    using parameter_type = potential_type::parameter_type;

    constexpr std::size_t N = 10000;
    constexpr real_type   h = 1e-6;
    constexpr real_type tol = 1e-6;

    const real_type sigma   = 3.0;
    const real_type epsilon = 1.0;
    const real_type omega   = 2.0;

    potential_type potential;

    const real_type x_min = 1.122462 * sigma;
    const real_type x_max = 1.122462 * sigma + omega;

    mjolnir::test::check_potential(potential,
            parameter_type{sigma, epsilon, omega, 1.0 / omega},
            x_min, x_max, tol, h, N);
}

BOOST_AUTO_TEST_CASE(iSoLF_float)
{
    using real_type = float;
    using potential_type = mjolnir::iSoLFAttractivePotential<real_type>;
    using parameter_type = potential_type::parameter_type;

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 0.002f;
    constexpr real_type tol = 0.004f;

    const real_type sigma   = 3.0f;
    const real_type epsilon = 1.0f;
    const real_type omega   = 2.0f;

    potential_type potential;

    const real_type x_min = 1.122462f * sigma;
    const real_type x_max = 1.122462f * sigma + omega;

    mjolnir::test::check_potential(potential,
            parameter_type{sigma, epsilon, omega, 1.0f / omega},
            x_min, x_max, tol, h, N);
}
