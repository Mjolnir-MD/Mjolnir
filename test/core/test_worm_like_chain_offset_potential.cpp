#define BOOST_TEST_MODULE "test_worm_like_chain_offset_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/check_potential.hpp>
#include <mjolnir/forcefield/local/WormLikeChainOffsetPotential.hpp>

BOOST_AUTO_TEST_CASE(WormLikeChainOffset_double)
{
    mjolnir::LoggerManager::set_default_logger("test_worm_like_chain_offset_potential.log");

    using real_type = double;
    constexpr std::size_t N = 900;
    constexpr real_type   h = 1e-6;
    constexpr real_type tol = 1e-6;
    const real_type       p = 3.9;
    const real_type      lc = 19.0;
    const real_type      l0 = 5.0;

    mjolnir::WormLikeChainOffsetPotential<real_type> pot(p, lc, l0);

    const real_type x_min = l0 + lc * 0.1;
    const real_type x_max = l0 + lc * 0.9;

    mjolnir::test::check_potential(pot, x_min, x_max, tol, h, N);
}

BOOST_AUTO_TEST_CASE(WormLikeChainOffset_float)
{
    using real_type = float;
    constexpr std::size_t N = 900;
    constexpr real_type   h = 1e-3;
    constexpr real_type tol = 1e-2;
    const real_type       p = 3.9;
    const real_type      lc = 19.0;
    const real_type      l0 = 5.0;

    mjolnir::WormLikeChainOffsetPotential<real_type> pot(p, lc, l0);

    const real_type x_min = l0 + lc * 0.1;
    const real_type x_max = l0 + lc * 0.9;

    mjolnir::test::check_potential(pot, x_min, x_max, tol, h, N);
}
