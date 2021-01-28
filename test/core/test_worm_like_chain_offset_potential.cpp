#define BOOST_TEST_MODULE "test_worm_like_chain_offset_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/forcefield/local/WormLikeChainOffsetPotential.hpp>

BOOST_AUTO_TEST_CASE(WormLikeChainOffset_double)
{
    mjolnir::LoggerManager::set_default_logger("test_worm_like_chain_offset_potential.log");

    using real_type = double;
    constexpr std::size_t N = 900;
    constexpr real_type   h = 1e-6;
    const real_type       p = 3.9;
    const real_type      lc = 19.0;
    const real_type      l0 = 5.0;

    mjolnir::WormLikeChainOffsetPotential<real_type> wormlikechainoffset(p, lc, l0);

    const real_type x_min = 0.0;
    const real_type x_max = (lc + l0) * 0.9;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=1; i<N; ++i)
    {
        const real_type x = x_min + dx * i;
        const real_type pot1 = wormlikechainoffset.potential(x + h);
        const real_type pot2 = wormlikechainoffset.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2.0 * h);
        const real_type deri = wormlikechainoffset.derivative(x);

        if(std::abs(dpot) > h && std::abs(deri) > h)
        {
            BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
        }
    }
}

BOOST_AUTO_TEST_CASE(WormLikeChainOffset_float)
{
    using real_type = float;
    constexpr std::size_t N = 900;
    constexpr real_type   h = 1e-3;
    const real_type       p = 3.9;
    const real_type      lc = 19.0;
    const real_type      l0 = 5.0;

    mjolnir::WormLikeChainOffsetPotential<real_type> wormlikechainoffset(p, lc, l0);

    const real_type x_min = 0.0;
    const real_type x_max = (lc + l0) * 0.9;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=1; i<N; ++i)
    {
        const real_type x = x_min + dx * i;
        const real_type pot1 = wormlikechainoffset.potential(x + h);
        const real_type pot2 = wormlikechainoffset.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2.0 * h);
        const real_type deri = wormlikechainoffset.derivative(x);

        if(std::abs(dpot) > h && std::abs(deri) > h)
        {
            BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
        }
    }
}
