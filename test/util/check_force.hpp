#ifndef MJOLNIR_TEST_UTIL_CHECK_FORCE_HPP
#define MJOLNIR_TEST_UTIL_CHECK_FORCE_HPP
#include <mjolnir/math/Vector.hpp>
#include <mjolnir/core/System.hpp>

#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/core/GlobalInteractionBase.hpp>
#include <type_traits>
#include <test/util/clear_system.hpp>

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

namespace mjolnir
{
namespace test
{

// Check if the sum of the forces is zero, only if the interaction is internal.

template<typename traitsT, typename Interaction>
void check_net_force(System<traitsT> sys,
                     const Interaction& interaction,
                     const typename traitsT::real_type tol)
{
    using real_type       = typename traitsT::real_type;
    using coordinate_type = typename traitsT::coordinate_type;
    clear_force(sys);

    sys.preprocess_forces();
    interaction.calc_force(sys);
    sys.postprocess_forces();

    coordinate_type f_tot = math::make_coordinate<coordinate_type>(0,0,0);
    for(const auto& f : sys.forces())
    {
        f_tot += f;
    }

    BOOST_TEST(mjolnir::math::X(f_tot) == real_type(0), boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Y(f_tot) == real_type(0), boost::test_tools::tolerance(tol));
    BOOST_TEST(mjolnir::math::Z(f_tot) == real_type(0), boost::test_tools::tolerance(tol));
    return;
}

// This checks force applied to each particle and the numerical difference of
// the corresponding energy

template<typename traitsT, typename Interaction>
void check_force(const System<traitsT>& init,
                 const Interaction& interaction,
                 const typename traitsT::real_type tol,
                 const typename traitsT::real_type dr,
                 const bool zero_net_force = true) // sum of external force (e.g. wall potential) is non-zero
{
    using real_type = typename traitsT::real_type;

    for(const auto& f : init.forces())
    {
        BOOST_TEST_REQUIRE(mjolnir::math::X(f) == real_type(0));
        BOOST_TEST_REQUIRE(mjolnir::math::Y(f) == real_type(0));
        BOOST_TEST_REQUIRE(mjolnir::math::Z(f) == real_type(0));
    }
    for(std::size_t i=0; i<9; ++i)
    {
        BOOST_TEST_REQUIRE(init.virial()[i] == real_type(0));
    }

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
            sys.preprocess_forces();
            interaction.calc_force(sys);
            sys.postprocess_forces();

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
            sys.preprocess_forces();
            interaction.calc_force(sys);
            sys.postprocess_forces();

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
            sys.preprocess_forces();
            interaction.calc_force(sys);
            sys.postprocess_forces();

            mjolnir::math::Z(sys.position(idx)) += dr;

            // calc U(x+dx)
            const auto E1 = interaction.calc_energy(sys);

            // central difference
            const auto dE = (E1 - E0) * 0.5;

            BOOST_TEST(-dE == dr * mjolnir::math::Z(sys.force(idx)),
                       boost::test_tools::tolerance(tol));
        }
    }

    // check if virial is not calculated in calc_force
    for(std::size_t i=0; i<9; ++i)
    {
        BOOST_TEST_REQUIRE(init.virial()[i] == real_type(0));
    }

    if(zero_net_force)
    {
        check_net_force(sys, interaction, tol);
    }
}

} // test
} // mjolnir
#endif// MJOLNIR_TEST_UTIL_CHECK_FORCE_HPP
