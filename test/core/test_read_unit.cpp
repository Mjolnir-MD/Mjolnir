#define BOOST_TEST_MODULE "test_unit_conversion"
#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <boost/mpl/list.hpp>
#include <test/util/make_empty_input.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/input/read_units.hpp>
#include <mjolnir/util/string.hpp>

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

    BOOST_TEST(Real(1.0) * cs::nm_to_angstrom() == Real(10.0), test::tolerance<Real>());
    BOOST_TEST(Real(1.0) * cs::angstrom_to_nm() == Real( 0.1), test::tolerance<Real>());
    BOOST_TEST(cs::nm_to_angstrom() * cs::angstrom_to_nm() == Real(1.0),
               test::tolerance<Real>());

    BOOST_TEST(Real(1.0)     * cs::nm_to_m() == Real(1e-9), test::tolerance<Real>());
    BOOST_TEST(Real(1.0)     * cs::m_to_nm() == Real(1e+9), test::tolerance<Real>());
    BOOST_TEST(cs::nm_to_m() * cs::m_to_nm() == Real(1.0),  test::tolerance<Real>());

    BOOST_TEST(Real(1.0) * cs::angstrom_to_m() == Real(1e-10),       test::tolerance<Real>());
    BOOST_TEST(Real(1.0) * cs::m_to_angstrom() == Real(1e+10),       test::tolerance<Real>());
    BOOST_TEST(cs::angstrom_to_m() * cs::m_to_angstrom() == Real(1.0),
               test::tolerance<Real>());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(unit_conversion_energy, Real, test_targets)
{
    using cs = mjolnir::unit::constants<Real>;

    BOOST_TEST(Real(1.0)      * cs::cal_to_J() == Real(4.184), test::tolerance<Real>());
    BOOST_TEST(Real(4.184)    * cs::J_to_cal() == Real(1.0),   test::tolerance<Real>());
    BOOST_TEST(cs::cal_to_J() * cs::J_to_cal() == Real(1.0),   test::tolerance<Real>());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_input_angstrom_kcalmol, Real, test_targets)
{
    mjolnir::LoggerManager::set_default_logger("test_unit.log");

    using namespace mjolnir::literals::string_literals;
    using phys = mjolnir::physics::constants<Real>;
    using unit = mjolnir::unit::constants<Real>;

    auto input = mjolnir::test::make_empty_input();
    auto& units = toml::get<toml::table>(input["units"]);

    units["length"] = "angstrom"_s;
    units["energy"] = "kcal/mol"_s;

    // update units and physical constants here
    phys::reset();
    mjolnir::read_units<mjolnir::SimulatorTraits<Real, mjolnir::UnlimitedBoundary>>(input);

    BOOST_TEST(phys::m_to_length() == unit::m_to_angstrom(), test::tolerance<Real>());
    BOOST_TEST(phys::length_to_m() == unit::angstrom_to_m(), test::tolerance<Real>());

    BOOST_TEST(phys::L_to_volume() == Real(1e-3) * // 1[L] = 1e-3 [m^3]
               unit::m_to_angstrom() * unit::m_to_angstrom() * unit::m_to_angstrom(),
               test::tolerance<Real>());
    BOOST_TEST(phys::volume_to_L() == Real(1e+3) * // 1[m^3] = 1e+3 [L]
               unit::angstrom_to_m() * unit::angstrom_to_m() * unit::angstrom_to_m(),
               test::tolerance<Real>());

    // NA [1/mol], kB [J/K], epsilon_0 [C^2/Jm]
    BOOST_TEST(phys::NA()   == unit::avogadro_constant(), test::tolerance<Real>());
    BOOST_TEST(phys::kB()   == unit::boltzmann_constant() * unit::J_to_cal() * Real(1e-3) *
                               phys::NA(), test::tolerance<Real>());
    BOOST_TEST(phys::eps0() == unit::vacuum_permittivity() / unit::elementary_charge() /
                               unit::elementary_charge() / unit::J_to_cal() / Real(1e-3) /
                               phys::NA() / unit::m_to_angstrom(),
                               test::tolerance<Real>());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_input_angstrom_kJmol, Real, test_targets)
{
    mjolnir::LoggerManager::set_default_logger("test_unit.log");

    using namespace mjolnir::literals::string_literals;
    using phys = mjolnir::physics::constants<Real>;
    using unit = mjolnir::unit::constants<Real>;


    auto input  = mjolnir::test::make_empty_input();
    auto& units = toml::get<toml::table>(input["units"]);

    units["length"] = "angstrom"_s;
    units["energy"] = "kJ/mol"_s;

    // reset units and physical constants here
    phys::reset();
    mjolnir::read_units<mjolnir::SimulatorTraits<Real, mjolnir::UnlimitedBoundary>>(input);

    BOOST_TEST(phys::m_to_length() == unit::m_to_angstrom(), test::tolerance<Real>());
    BOOST_TEST(phys::length_to_m() == unit::angstrom_to_m(), test::tolerance<Real>());

    BOOST_TEST(phys::L_to_volume() == Real(1e-3) * // 1[L] = 1e-3 [m^3]
               unit::m_to_angstrom() * unit::m_to_angstrom() * unit::m_to_angstrom(),
               test::tolerance<Real>());
    BOOST_TEST(phys::volume_to_L() == Real(1e+3) * // 1[m^3] = 1e+3 [L]
               unit::angstrom_to_m() * unit::angstrom_to_m() * unit::angstrom_to_m(),
               test::tolerance<Real>());

    // NA [1/mol], kB [J/K], epsilon_0 [C^2/Jm]
    BOOST_TEST(phys::NA()   == unit::avogadro_constant(), test::tolerance<Real>());
    BOOST_TEST(phys::kB()   == unit::boltzmann_constant() * Real(1e-3) *
                               phys::NA(), test::tolerance<Real>());
    BOOST_TEST(phys::eps0() == unit::vacuum_permittivity() / unit::elementary_charge() /
                               unit::elementary_charge()  / Real(1e-3) /
                               phys::NA() / unit::m_to_angstrom(),
                               test::tolerance<Real>());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_input_nm_kcalmol, Real, test_targets)
{
    mjolnir::LoggerManager::set_default_logger("test_unit.log");

    using namespace mjolnir::literals::string_literals;
    using phys = mjolnir::physics::constants<Real>;
    using unit = mjolnir::unit::constants<Real>;

    auto input = mjolnir::test::make_empty_input();
    auto& units = toml::get<toml::table>(input["units"]);

    units["length"] = "nm"_s;
    units["energy"] = "kcal/mol"_s;

    // update units and physical constants here
    phys::reset();
    mjolnir::read_units<mjolnir::SimulatorTraits<Real, mjolnir::UnlimitedBoundary>>(input);

    BOOST_TEST(phys::m_to_length() == unit::m_to_nm(), test::tolerance<Real>());
    BOOST_TEST(phys::length_to_m() == unit::nm_to_m(), test::tolerance<Real>());

    BOOST_TEST(phys::L_to_volume() == Real(1e-3) * // 1[L] = 1e-3 [m^3]
               unit::m_to_nm() * unit::m_to_nm() * unit::m_to_nm(),
               test::tolerance<Real>());
    BOOST_TEST(phys::volume_to_L() == Real(1e+3) * // 1[m^3] = 1e+3 [L]
               unit::nm_to_m() * unit::nm_to_m() * unit::nm_to_m(),
               test::tolerance<Real>());

    // NA [1/mol], kB [J/K], epsilon_0 [C^2/Jm]
    BOOST_TEST(phys::NA()   == unit::avogadro_constant(), test::tolerance<Real>());
    BOOST_TEST(phys::kB()   == unit::boltzmann_constant() * unit::J_to_cal() * Real(1e-3) *
                               phys::NA(), test::tolerance<Real>());
    BOOST_TEST(phys::eps0() == unit::vacuum_permittivity() / unit::elementary_charge() /
                               unit::elementary_charge() / unit::J_to_cal() / Real(1e-3) /
                               phys::NA() / unit::m_to_nm(),
                               test::tolerance<Real>());
}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_input_nm_kJmol, Real, test_targets)
{
    mjolnir::LoggerManager::set_default_logger("test_unit.log");

    using namespace mjolnir::literals::string_literals;
    using phys = mjolnir::physics::constants<Real>;
    using unit = mjolnir::unit::constants<Real>;

    auto input  = mjolnir::test::make_empty_input();
    auto& units = toml::get<toml::table>(input["units"]);

    units["length"] = "nm"_s;
    units["energy"] = "kJ/mol"_s;

    // update units and physical constants here
    phys::reset();
    mjolnir::read_units<mjolnir::SimulatorTraits<Real, mjolnir::UnlimitedBoundary>>(input);

    BOOST_TEST(phys::m_to_length() == unit::m_to_nm(), test::tolerance<Real>());
    BOOST_TEST(phys::length_to_m() == unit::nm_to_m(), test::tolerance<Real>());

    BOOST_TEST(phys::L_to_volume() == Real(1e-3) * // 1[L] = 1e-3 [m^3]
               unit::m_to_nm() * unit::m_to_nm() * unit::m_to_nm(),
               test::tolerance<Real>());
    BOOST_TEST(phys::volume_to_L() == Real(1e+3) * // 1[m^3] = 1e+3 [L]
               unit::nm_to_m() * unit::nm_to_m() * unit::nm_to_m(),
               test::tolerance<Real>());

    // NA [1/mol], kB [J/K], epsilon_0 [C^2/Jm]
    BOOST_TEST(phys::NA()   == unit::avogadro_constant(), test::tolerance<Real>());
    BOOST_TEST(phys::kB()   == unit::boltzmann_constant() * Real(1e-3) *
                               phys::NA(), test::tolerance<Real>());
    BOOST_TEST(phys::eps0() == unit::vacuum_permittivity() / unit::elementary_charge() /
                               unit::elementary_charge()  / Real(1e-3) /
                               phys::NA() / unit::m_to_nm(),
                               test::tolerance<Real>());
}
