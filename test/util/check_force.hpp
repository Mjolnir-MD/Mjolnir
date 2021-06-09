#ifndef MJOLNIR_TEST_UTIL_CHECK_FORCE_HPP
#define MJOLNIR_TEST_UTIL_CHECK_FORCE_HPP
#include <mjolnir/math/Vector.hpp>
#include <mjolnir/core/System.hpp>

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

template<typename traitsT, typename Interaction>
void check_force(const System<traitsT>& init,
                 const Interaction& interaction,
                 const typename traitsT::real_type tol,
                 const typename traitsT::real_type dr)
{
    System<traitsT> sys(init);

    for(std::size_t idx=0; idx<sys.size(); ++idx)
    {
        {
            // ----------------------------------------------------------------
            // reset positions
            sys = init;

            // calc U(x-dx)
            const auto E0 = interaction.calc_energy(sys);

            mjolnir::math::X(sys.position(idx)) += dr;

            // calc F(x)
            interaction.calc_force(sys);

            mjolnir::math::X(sys.position(idx)) += dr;

            // calc U(x+dx)
            const auto E1 = interaction.calc_energy(sys);

            // central difference
            const auto dE = (E1 - E0) * 0.5;

            BOOST_TEST(-dE == dr * mjolnir::math::X(sys.force(idx)),
                       boost::test_tools::tolerance(tol));
        }
        {
            // ----------------------------------------------------------------
            // reset positions
            sys = init;

            // calc U(x-dx)
            const auto E0 = interaction.calc_energy(sys);

            mjolnir::math::Y(sys.position(idx)) += dr;

            // calc F(x)
            interaction.calc_force(sys);

            mjolnir::math::Y(sys.position(idx)) += dr;

            // calc U(x+dx)
            const auto E1 = interaction.calc_energy(sys);

            // central difference
            const auto dE = (E1 - E0) * 0.5;

            BOOST_TEST(-dE == dr * mjolnir::math::Y(sys.force(idx)),
                       boost::test_tools::tolerance(tol));
        }
        {
            // ----------------------------------------------------------------
            // reset positions
            sys = init;

            // calc U(x-dx)
            const auto E0 = interaction.calc_energy(sys);

            mjolnir::math::Z(sys.position(idx)) += dr;

            // calc F(x)
            interaction.calc_force(sys);

            mjolnir::math::Z(sys.position(idx)) += dr;

            // calc U(x+dx)
            const auto E1 = interaction.calc_energy(sys);

            // central difference
            const auto dE = (E1 - E0) * 0.5;

            BOOST_TEST(-dE == dr * mjolnir::math::Z(sys.force(idx)),
                       boost::test_tools::tolerance(tol));
        }
    }
}

} // test
} // mjolnir
#endif// MJOLNIR_TEST_UTIL_CHECK_FORCE_HPP
