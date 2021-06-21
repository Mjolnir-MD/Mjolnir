#define BOOST_TEST_MODULE "test_wca_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/check_potential.hpp>
#include <mjolnir/forcefield/global/WCAPotential.hpp>

BOOST_AUTO_TEST_CASE(WCA_double)
{
    using real_type = double;
    using potential_type = mjolnir::WCAPotential<real_type>;
    using parameter_type = potential_type::parameter_type;

    constexpr static std::size_t N = 10000;
    constexpr static real_type   h = 1e-6;
    constexpr static real_type tol = 1e-6;

    const real_type sigma   = 3.0;
    const real_type epsilon = 1.0;

    potential_type potential;

    const real_type x_min = 0.8 * sigma;
    const real_type x_max = potential.cutoff_ratio() * sigma - 2 * h;

    mjolnir::test::check_potential(potential, parameter_type{sigma, epsilon},
                                   x_min, x_max, tol, h, N);
}

BOOST_AUTO_TEST_CASE(WCA_float)
{
    using real_type = double;
    using potential_type = mjolnir::WCAPotential<real_type>;
    using parameter_type = potential_type::parameter_type;

    constexpr static std::size_t N = 1000;
    constexpr static real_type   h = 0.002f;
    constexpr static real_type tol = 0.005f;

    const real_type sigma   = 3.0f;
    const real_type epsilon = 1.0f;

    potential_type potential;

    const real_type x_min = 0.8f * sigma;
    const real_type x_max = potential.cutoff_ratio() * sigma - 2 * h;

    mjolnir::test::check_potential(potential, parameter_type{sigma, epsilon},
                                   x_min, x_max, tol, h, N);
}


