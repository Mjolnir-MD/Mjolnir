#define BOOST_TEST_MODULE "test_throw_exception"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/util/throw_exception.hpp>
#include <stdexcept>

BOOST_AUTO_TEST_CASE(test_what_message)
{
    {
        std::string what;
        try
        {
            mjolnir::throw_exception<std::runtime_error>("st", "ring");
        }
        catch(std::runtime_error const& re)
        {
            what = re.what();
        }

        std::ostringstream oss;
        oss << "st" << "ring";

        BOOST_CHECK_EQUAL(oss.str(), what);
    }
    {
        const int i = 42;
        std::string what;
        try
        {
            mjolnir::throw_exception<std::runtime_error>("integer i = ", i);
        }
        catch(std::runtime_error const& re)
        {
            what = re.what();
        }

        std::ostringstream oss;
        oss << "integer i = " << i;

        BOOST_CHECK_EQUAL(oss.str(), what);
    }
    {
        const int    i = 42;
        const double d = 3.14;
        std::string what;
        try
        {
            mjolnir::throw_exception<std::runtime_error>(
                    "integer i = ", i, ", double d = ", d);
        }
        catch(std::runtime_error const& re)
        {
            what = re.what();
        }

        std::ostringstream oss;
        oss << "integer i = " << i << ", double d = " << d;

        BOOST_CHECK_EQUAL(oss.str(), what);
    }
}

BOOST_AUTO_TEST_CASE(test_empty_message)
{
    {
        std::string what;
        try
        {
            mjolnir::throw_exception<std::runtime_error>();
        }
        catch(std::runtime_error const& re)
        {
            what = re.what();
        }

        BOOST_CHECK(what.empty());
    }
}
