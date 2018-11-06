#define BOOST_TEST_MODULE "test_unit_conversion"
#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/input/read_units.hpp>

typedef boost::mpl::list<double, float> test_targets;

namespace test
{
template<typename T>
decltype(boost::test_tools::tolerance(std::declval<T>())) tolerance();
template<>
decltype(boost::test_tools::tolerance(std::declval<double>()))
tolerance<double>() {return boost::test_tools::tolerance(1e-8);}
template<>
decltype(boost::test_tools::tolerance(std::declval<float>()))
tolerance<float>()  {return boost::test_tools::tolerance(1e-4f);}
} // test

BOOST_AUTO_TEST_CASE_TEMPLATE(unit_conversion_length, Real, test_targets)
{
    using cs = mjolnir::unit::constants<Real>;

    BOOST_TEST(Real(1.0) * cs::nm_to_angstrom == Real(10.0), test::tolerance<Real>());
    BOOST_TEST(Real(1.0) * cs::angstrom_to_nm == Real( 0.1), test::tolerance<Real>());
    BOOST_TEST(cs::nm_to_angstrom * cs::angstrom_to_nm == Real(1.0),
               test::tolerance<Real>());

    BOOST_TEST(Real(1.0)   * cs::nm_to_m == Real(1e-9), test::tolerance<Real>());
    BOOST_TEST(Real(1.0)   * cs::m_to_nm == Real(1e+9), test::tolerance<Real>());
    BOOST_TEST(cs::nm_to_m * cs::m_to_nm == Real(1.0),  test::tolerance<Real>());

    BOOST_TEST(Real(1.0) * cs::angstrom_to_m == Real(1e-10), test::tolerance<Real>());
    BOOST_TEST(Real(1.0) * cs::m_to_angstrom == Real(1e+10), test::tolerance<Real>());
    BOOST_TEST(cs::angstrom_to_m * cs::m_to_angstrom == Real(1.0), test::tolerance<Real>());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(unit_conversion_energy, Real, test_targets)
{
    using cs = mjolnir::unit::constants<Real>;

    BOOST_TEST(Real(1.0)    * cs::cal_to_J == Real(4.1868), test::tolerance<Real>());
    BOOST_TEST(Real(4.1868) * cs::J_to_cal == Real(1.0), test::tolerance<Real>());
    BOOST_TEST(cs::cal_to_J * cs::J_to_cal == Real(1.0), test::tolerance<Real>());
}
