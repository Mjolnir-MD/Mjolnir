#define BOOST_TEST_MODULE "test_gaussian_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/check_potential.hpp>
#include <mjolnir/forcefield/local/GaussianPotential.hpp>

BOOST_AUTO_TEST_CASE(Gaussian_double)
{
    using real_type = double;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 1e-6;
    constexpr real_type tol = 1e-6;
    const real_type e  = 2.0;
    const real_type w  = 0.15;
    const real_type r0 = 7.0;

    mjolnir::GaussianPotential<real_type> gaussian(e, w, r0);

    const real_type x_min = r0 - 5.0 * w;
    const real_type x_max = r0 + 5.0 * w;

    mjolnir::test::check_potential(gaussian, x_min, x_max, tol, h, N);
}

BOOST_AUTO_TEST_CASE(Gaussian_float)
{
    using real_type = float;
    constexpr std::size_t N = 100;
    constexpr real_type   h = 1e-3f;
    constexpr real_type tol = 1e-3f;
    const real_type e  = 2.0;
    const real_type w  = 0.15;
    const real_type r0 = 7.0;

    mjolnir::GaussianPotential<real_type> gaussian(e, w, r0);

    const real_type x_min = r0 - 5.0f * w;
    const real_type x_max = r0 + 5.0f * w;
    mjolnir::test::check_potential(gaussian, x_min, x_max, tol, h, N);
}

BOOST_AUTO_TEST_CASE(Gaussian_cutoff_double)
{
    using real_type = double;
    constexpr real_type   h = 1e-6;

    {
        const real_type e  = 2.0;
        const real_type w  = 0.15;
        const real_type r0 = 7.0;

        mjolnir::GaussianPotential<real_type> gaussian(e, w, r0);

        const auto rc = gaussian.cutoff();

        BOOST_TEST(std::isfinite(rc));
        BOOST_TEST(std::abs(gaussian.potential(rc)) ==
                   mjolnir::math::abs_tolerance<real_type>(),
                   boost::test_tools::tolerance(h));
    }
    {
        const real_type e  = -2.0;
        const real_type w  = 0.15;
        const real_type r0 = 7.0;

        mjolnir::GaussianPotential<real_type> gaussian(e, w, r0);

        const auto rc = gaussian.cutoff();

        BOOST_TEST(std::isfinite(rc));
        BOOST_TEST(std::abs(gaussian.potential(rc)) ==
                   mjolnir::math::abs_tolerance<real_type>(),
                   boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(Gaussian_cutoff_float)
{
    using real_type = float;
    constexpr real_type   h = 1e-3;

    {
        const real_type e  = 2.0;
        const real_type w  = 0.15;
        const real_type r0 = 7.0;

        mjolnir::GaussianPotential<real_type> gaussian(e, w, r0);

        const auto rc = gaussian.cutoff();

        BOOST_TEST(std::isfinite(rc));
        BOOST_TEST(gaussian.potential(rc) ==
                   mjolnir::math::abs_tolerance<real_type>(),
                   boost::test_tools::tolerance(h));
    }
    {
        const real_type e  = -2.0;
        const real_type w  = 0.15;
        const real_type r0 = 7.0;

        mjolnir::GaussianPotential<real_type> gaussian(e, w, r0);

        const auto rc = gaussian.cutoff();

        BOOST_TEST(std::isfinite(rc));
        BOOST_TEST(std::abs(gaussian.potential(rc)) ==
                   mjolnir::math::abs_tolerance<real_type>(),
                   boost::test_tools::tolerance(h));
    }
}
