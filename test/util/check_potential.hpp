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
    using real_type = typename Potential::real_type;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=1; i<N; ++i)
    {
        const real_type x = x_min + i * dx;
        const real_type pot1 = pot.potential(x + h);
        const real_type pot2 = pot.potential(x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);

        const real_type deri = pot.derivative(x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(tol));
    }
}

template<typename Potential>
void check_potential(const Potential& pot,
        const typename Potential::parameter_type& para,
        const typename Potential::real_type x_min,
        const typename Potential::real_type x_max,
        const typename Potential::real_type tol,
        const typename Potential::real_type h,
        const std::size_t N)
{
    using real_type = typename Potential::real_type;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=1; i<N; ++i)
    {
        const real_type x = x_min + i * dx;
        const real_type pot1 = pot.potential(x + h, para);
        const real_type pot2 = pot.potential(x - h, para);
        const real_type dpot = (pot1 - pot2) / (2 * h);

        const real_type deri = pot.derivative(x, para);

        BOOST_TEST_MESSAGE("x = " << x << ", V(x+h) = " << pot1 << ", V(x-h) = " << pot2
                           << ", dV/h = " << dpot << ", dV/dx = " << deri);
        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(tol));
    }
}

} // test
} // mjolnir
#endif// MJOLNIR_TEST_UTIL_CHECK_POTENTIAL_HPP
