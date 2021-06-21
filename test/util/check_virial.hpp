#ifndef MJOLNIR_TEST_UTIL_CHECK_VIRIAL_HPP
#define MJOLNIR_TEST_UTIL_CHECK_VIRIAL_HPP
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

// To check the virial, it compares virial calculated from forcefield and r*f.
// This relationship holds only when the interacting pairs do not stick out of
// the box. Thus it checks the virial only if the system is under the unlimited
// boundary condition.
template<typename traitsT, typename Interaction>
void check_virial(System<traitsT>& sys,
                  const Interaction& interaction,
                  const typename traitsT::real_type tol)
{
    using coordinate_type = typename traitsT::coordinate_type;
    using matrix33_type = typename traitsT::matrix33_type;

    if(is_unlimited_boundary<traitsT>::value)
    {
        sys.virial() = matrix33_type(0,0,0, 0,0,0, 0,0,0);
        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            sys.force(idx) = math::make_coordinate<coordinate_type>(0,0,0);
        }
        sys.preprocess_forces();
        interaction.calc_force(sys);
        sys.postprocess_forces();

        matrix33_type vir(0,0,0, 0,0,0, 0,0,0);
        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            vir += math::tensor_product(sys.position(idx), sys.force(idx));
        }

        BOOST_TEST(sys.virial()(0,0) == vir(0,0), boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.virial()(0,1) == vir(0,1), boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.virial()(0,2) == vir(0,2), boost::test_tools::tolerance(tol));

        BOOST_TEST(sys.virial()(1,0) == vir(1,0), boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.virial()(1,1) == vir(1,1), boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.virial()(1,2) == vir(1,2), boost::test_tools::tolerance(tol));

        BOOST_TEST(sys.virial()(2,0) == vir(2,0), boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.virial()(2,1) == vir(2,1), boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.virial()(2,2) == vir(2,2), boost::test_tools::tolerance(tol));
    }
}

} // test
} // mjolnir
#endif// MJOLNIR_TEST_UTIL_CHECK_FORCE_HPP
