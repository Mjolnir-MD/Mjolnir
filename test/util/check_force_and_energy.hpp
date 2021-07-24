#ifndef MJOLNIR_TEST_UTIL_CHECK_FORCE_AND_ENERGY_HPP
#define MJOLNIR_TEST_UTIL_CHECK_FORCE_AND_ENERGY_HPP
#include <mjolnir/math/Vector.hpp>
#include <mjolnir/core/System.hpp>

#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/core/GlobalInteractionBase.hpp>
#include <type_traits>

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

namespace mjolnir
{
namespace test
{

// This checks the consistency between `calc_force` and `calc_force_and_energy`.

template<typename traitsT, typename Interaction>
void check_force_and_energy(System<traitsT> ref,
                            const Interaction& interaction,
                            const typename traitsT::real_type tol)
{
    using coordinate_type = typename traitsT::coordinate_type;
    using matrix33_type   = typename traitsT::matrix33_type;

    for(std::size_t i=0; i<ref.size(); ++i)
    {
        ref.force(i) = math::make_coordinate<coordinate_type>(0,0,0);
    }
    ref.virial() = matrix33_type(0,0,0, 0,0,0, 0,0,0);

    System<traitsT> sys(ref);

    ref.preprocess_forces();
    interaction.calc_force(ref);
    ref.postprocess_forces();
    const auto ref_ene = interaction.calc_energy(ref);

    sys.preprocess_forces();
    const auto ene = interaction.calc_force_and_energy(sys);
    sys.postprocess_forces();

    // ------------------------------------------------------------------------
    // energy, forces are the same. virial is not calculated

    BOOST_TEST(ref_ene == ene, boost::test_tools::tolerance(tol));

    for(std::size_t idx=0; idx<sys.size(); ++idx)
    {
        BOOST_TEST(math::X(sys.force(idx)) == math::X(ref.force(idx)), boost::test_tools::tolerance(tol));
        BOOST_TEST(math::Y(sys.force(idx)) == math::Y(ref.force(idx)), boost::test_tools::tolerance(tol));
        BOOST_TEST(math::Z(sys.force(idx)) == math::Z(ref.force(idx)), boost::test_tools::tolerance(tol));
    }
    for(std::size_t i=0; i<9; ++i)
    {
        BOOST_TEST(sys.virial()[i] == 0.0, boost::test_tools::tolerance(tol));
    }

    return;
}

} // test
} // mjolnir
#endif// MJOLNIR_TEST_UTIL_CHECK_FORCE_HPP
