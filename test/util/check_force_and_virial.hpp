#ifndef MJOLNIR_TEST_UTIL_CHECK_FORCE_AND_VIRIAL_HPP
#define MJOLNIR_TEST_UTIL_CHECK_FORCE_AND_VIRIAL_HPP
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

// This checks the consistency of forces between `calc_force` and
// `calc_force_and_virial`.

template<typename traitsT, typename Interaction>
void check_force_and_virial(System<traitsT> ref,
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

    sys.preprocess_forces();
    interaction.calc_force_and_virial(sys);
    sys.postprocess_forces();

    for(std::size_t idx=0; idx<sys.size(); ++idx)
    {
        BOOST_TEST(math::X(sys.force(idx)) == math::X(ref.force(idx)), boost::test_tools::tolerance(tol));
        BOOST_TEST(math::Y(sys.force(idx)) == math::Y(ref.force(idx)), boost::test_tools::tolerance(tol));
        BOOST_TEST(math::Z(sys.force(idx)) == math::Z(ref.force(idx)), boost::test_tools::tolerance(tol));
    }
}

} // test
} // mjolnir
#endif// MJOLNIR_TEST_UTIL_CHECK_FORCE_HPP
