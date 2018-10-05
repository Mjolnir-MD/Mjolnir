#define BOOST_TEST_MODULE "test_unit_conversion"
#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <mjolnir/core/Unit.hpp>

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

BOOST_AUTO_TEST_CASE_TEMPLATE(Unit_conversion_length, Real, test_targets)
{
    //XXX here you need to convert floating-point literal explicitly because
    //    it will be expanded to template code that may cast float into double.
    using namespace mjolnir;

    /* um -> um */{
        unit::quantity<Real, unit::micro_meter> um(Real(1));
        unit::quantity<Real, unit::micro_meter> um1(um);

        BOOST_TEST(um1.value == Real(1), test::tolerance<Real>());
    }
    /* um -> nm */{
        unit::quantity<Real, unit::micro_meter> um(Real(1));
        unit::quantity<Real, unit::nano_meter > nm(um);

        BOOST_TEST(nm.value == Real(1000), test::tolerance<Real>());
    }
    /* um -> angstrom */{
        unit::quantity<Real, unit::micro_meter> um(Real(1));
        unit::quantity<Real, unit::angstrom   > A(um);

        BOOST_TEST(A.value == Real(10000), test::tolerance<Real>());
    }
    /* um -> pm */{
        unit::quantity<Real, unit::micro_meter> um(Real(1));
        unit::quantity<Real, unit::pico_meter > pm(um);

        BOOST_TEST(pm.value == Real(1000000), test::tolerance<Real>());
    }


    /* nm -> um */ {
        unit::quantity<Real, unit::nano_meter>  nm(Real(1));
        unit::quantity<Real, unit::micro_meter> um(nm);

        BOOST_TEST(um.value == Real(0.001), test::tolerance<Real>());
    }
    /* nm -> nm */ {
        unit::quantity<Real, unit::nano_meter> nm(Real(1));
        unit::quantity<Real, unit::nano_meter> nm1(nm);

        BOOST_TEST(nm1.value == Real(1), test::tolerance<Real>());
    }
    /* nm -> angstrom */ {
        unit::quantity<Real, unit::nano_meter> nm(Real(1));
        unit::quantity<Real, unit::angstrom  > A(nm);

        BOOST_TEST(A.value == Real(10.0), test::tolerance<Real>());
    }
    /* nm -> pm */ {
        unit::quantity<Real, unit::nano_meter> nm(Real(1));
        unit::quantity<Real, unit::pico_meter> pm(nm);

        BOOST_TEST(pm.value == Real(1000.0), test::tolerance<Real>());
    }

    /* A -> um */ {
        unit::quantity<Real, unit::angstrom>    A(Real(1));
        unit::quantity<Real, unit::micro_meter> um(A);

        BOOST_TEST(um.value == Real(0.0001), test::tolerance<Real>());
    }
    /* A -> nm */ {
        unit::quantity<Real, unit::angstrom>   A(Real(1));
        unit::quantity<Real, unit::nano_meter> nm(A);

        BOOST_TEST(nm.value == Real(0.1), test::tolerance<Real>());
    }
    /* A -> angstrom */ {
        unit::quantity<Real, unit::angstrom> A(Real(1));
        unit::quantity<Real, unit::angstrom> A1(A);

        BOOST_TEST(A1.value == Real(1.0), test::tolerance<Real>());
    }
    /* A -> pm */ {
        unit::quantity<Real, unit::angstrom>   A(Real(1));
        unit::quantity<Real, unit::pico_meter> pm(A);

        BOOST_TEST(pm.value == Real(100.0), test::tolerance<Real>());
    }

    /* pm -> um */ {
        unit::quantity<Real, unit::pico_meter>   pm(Real(1));
        unit::quantity<Real, unit::micro_meter> um(pm);

        BOOST_TEST(um.value == Real(0.000001), test::tolerance<Real>());
    }
    /* pm -> nm */ {
        unit::quantity<Real, unit::pico_meter>   pm(Real(1));
        unit::quantity<Real, unit::nano_meter> nm(pm);

        BOOST_TEST(nm.value == Real(0.001), test::tolerance<Real>());
    }
    /* pm -> angstrom */ {
        unit::quantity<Real, unit::pico_meter> pm(Real(1));
        unit::quantity<Real, unit::angstrom> A(pm);

        BOOST_TEST(A.value == Real(0.01), test::tolerance<Real>());
    }
    /* pm -> pm */ {
        unit::quantity<Real, unit::pico_meter>   pm(Real(1));
        unit::quantity<Real, unit::pico_meter> pm1(pm);

        BOOST_TEST(pm.value == Real(1.0), test::tolerance<Real>());
    }
}
