#define BOOST_TEST_MODULE "test_hard_core_excluded_volume_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>
#include <mjolnir/forcefield/global/HardCoreExcludedVolumePotential.hpp>

BOOST_AUTO_TEST_CASE(HCEXV_double)
{
    using real_type = double;
    using potential_type = mjolnir::HardCoreExcludedVolumePotential<real_type>;
    using parameter_type = potential_type::parameter_type;

    constexpr static std::size_t N = 10000;
    constexpr static real_type   h = 1e-6;
    constexpr static real_type tol = 1e-6;

    const real_type epsilon = 1.0;
    const real_type hard_core_radius    = 10.0;
    const real_type soft_shell_thickness = 3.0;

    potential_type potential(2.0, epsilon);

    const parameter_type params{soft_shell_thickness, hard_core_radius};

    const real_type x_min = 0.8 * soft_shell_thickness + hard_core_radius;
    const real_type x_max = potential.absolute_cutoff(params);

    mjolnir::test::check_potential(potential, params,
                                   x_min, x_max, tol, h, N);
}

BOOST_AUTO_TEST_CASE(HCEXV_float)
{
    using real_type = float;
    using potential_type = mjolnir::HardCoreExcludedVolumePotential<real_type>;
    using parameter_type = potential_type::parameter_type;

    constexpr static std::size_t N = 1000;
    constexpr static real_type   h = 0.001f;
    constexpr static real_type tol = 0.002f;

    const real_type epsilon = 1.0f;
    const real_type hard_core_radius    = 10.0f;
    const real_type soft_shell_thickness = 3.0f;

    potential_type potential(2.0f, epsilon);

    const parameter_type params{soft_shell_thickness, hard_core_radius};

    const real_type x_min = 0.8f * soft_shell_thickness + hard_core_radius;
    const real_type x_max = potential.absolute_cutoff(params);

    mjolnir::test::check_potential(potential, params,
                                   x_min, x_max, tol, h, N);
}
