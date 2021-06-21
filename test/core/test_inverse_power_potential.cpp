#define BOOST_TEST_MODULE "test_inverse_power_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/forcefield/global/InversePowerPotential.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

BOOST_AUTO_TEST_CASE(InversePower_double)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_inverse_power_potential.log");

    using real_type = double;
    using potential_type = mjolnir::InversePowerPotential<real_type>;
    using parameter_type = typename potential_type::parameter_type;

    constexpr std::size_t N = 10000;
    constexpr real_type   h = 1e-6;
    constexpr real_type tol = 1e-6;

    const real_type sigma   = 3.0;
    const real_type epsilon = 1.0;
    const std::int32_t n    = 5;

    potential_type potential(potential_type::default_cutoff(n), epsilon, n);

    const real_type x_min = 0.8 * sigma;
    const real_type x_max = potential.absolute_cutoff(parameter_type{sigma});

    mjolnir::test::check_potential(potential, parameter_type{sigma},
                                   x_min, x_max, tol, h, N);
}

BOOST_AUTO_TEST_CASE(InversePower_float)
{
    mjolnir::LoggerManager::set_default_logger(
            "test_inverse_power_potential.log");

    using real_type = float;
    using potential_type = mjolnir::InversePowerPotential<real_type>;
    using parameter_type = typename potential_type::parameter_type;

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 0.001f;
    constexpr real_type tol = 0.002f;

    const real_type sigma   = 3.0f;
    const real_type epsilon = 1.0f;
    const std::int32_t n    = 5;

    potential_type potential(potential_type::default_cutoff(n), epsilon, n);

    const real_type x_min = 0.8f * sigma;
    const real_type x_max = potential.absolute_cutoff(parameter_type{sigma});

    mjolnir::test::check_potential(potential, parameter_type{sigma},
                                   x_min, x_max, tol, h, N);
}
