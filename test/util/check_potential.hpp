#ifndef MJOLNIR_TEST_UTIL_CHECK_POTENTIAL_HPP
#define MJOLNIR_TEST_UTIL_CHECK_POTENTIAL_HPP

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

namespace mjolnir
{
namespace test
{

// This checks force applied to each particle and the numerical difference of
// the corresponding energy

template<typename Potential>
void check_potential(const Potential& pot,
        const typename Potential::real_type x_min,
        const typename Potential::real_type x_max,
        const typename Potential::real_type tol,
        const typename Potential::real_type h,
        const std::size_t N)
{
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x = x_min + i * dx;
        const real_type pot1 = pot.potential(x + h);
        const real_type pot2 = pot.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);

        const real_type deri = pot.derivative(x);

        if(std::abs(pot1 / pot2 - 1.0) < h)
        {
            // pot1 and pot2 are almost the same, thus dpot ~ 0.0.
            // to avoid numerical error in dpot, here it checks `deri ~ 0.0`.
            BOOST_TEST(deri == 0.0, boost::test_tools::tolerance(h));
        }
        else
        {
            BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
        }
    }
}

} // test
} // mjolnir
#endif// MJOLNIR_TEST_UTIL_CHECK_POTENTIAL_HPP
