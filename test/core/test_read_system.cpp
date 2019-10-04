#define BOOST_TEST_MODULE "test_read_system"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/input/read_system.hpp>

BOOST_AUTO_TEST_CASE(read_system_with_unlimited_boundary)
{
    mjolnir::LoggerManager::set_default_logger("test_read_system.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    constexpr real_type tol = 1e-8;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            [[systems]]
            boundary_shape = {}
            attributes = {test_attr = 3.14}
            particles  = [
                {m = 1.0, pos = [1.0, 2.0, 3.0], vel = [1.1, 2.1, 3.1], name = "testnm", group = "testgrp"},
                {m = 2.0, pos = [1.2, 2.2, 3.2], vel = [1.3, 2.3, 3.3], name = "testnm", group = "testgrp"},
                {m = 3.0, pos = [1.4, 2.4, 3.4], vel = [1.5, 2.5, 3.5]}
            ]
        )"_toml;

        const auto sys = mjolnir::read_system<traits_type>(v, 0);
        BOOST_TEST(sys.size() == 3u);

        BOOST_TEST(sys.mass(0) == 1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.mass(1) == 2.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.mass(2) == 3.0, boost::test_tools::tolerance(tol));

        BOOST_TEST(sys.rmass(0) == 1.0 / 1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.rmass(1) == 1.0 / 2.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.rmass(2) == 1.0 / 3.0, boost::test_tools::tolerance(tol));

        BOOST_TEST(sys.position(0).at(0) == 1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.position(0).at(1) == 2.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.position(0).at(2) == 3.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.position(1).at(0) == 1.2, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.position(1).at(1) == 2.2, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.position(1).at(2) == 3.2, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.position(2).at(0) == 1.4, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.position(2).at(1) == 2.4, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.position(2).at(2) == 3.4, boost::test_tools::tolerance(tol));

        BOOST_TEST(sys.velocity(0).at(0) == 1.1, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.velocity(0).at(1) == 2.1, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.velocity(0).at(2) == 3.1, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.velocity(1).at(0) == 1.3, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.velocity(1).at(1) == 2.3, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.velocity(1).at(2) == 3.3, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.velocity(2).at(0) == 1.5, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.velocity(2).at(1) == 2.5, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.velocity(2).at(2) == 3.5, boost::test_tools::tolerance(tol));

        BOOST_TEST(sys.force(0).at(0) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(0).at(1) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(0).at(2) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(1).at(0) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(1).at(1) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(1).at(2) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(2).at(0) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(2).at(1) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(2).at(2) == 0.0, boost::test_tools::tolerance(tol));

        BOOST_TEST(sys.name(0) == "testnm");
        BOOST_TEST(sys.name(1) == "testnm");
        BOOST_TEST(sys.name(2) == "X");

        BOOST_TEST(sys.group(0) == "testgrp");
        BOOST_TEST(sys.group(1) == "testgrp");
        BOOST_TEST(sys.group(2) == "NONE");

        BOOST_TEST(sys.has_attribute("test_attr"));
        BOOST_TEST(sys.attribute("test_attr") == 3.14, boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_system_with_cuboidal_periodic_boundary)
{
    mjolnir::LoggerManager::set_default_logger("test_read_system.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::CuboidalPeriodicBoundary>;
    constexpr real_type tol = 1e-8;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            [[systems]]
            boundary_shape.upper = [10.0, 11.0, 12.0]
            boundary_shape.lower = [ 0.0,  1.0,  2.0]
            attributes = {test_attr = 3.14}
            particles  = [
                {m = 1.0, pos = [1.0, 2.0, 3.0], vel = [1.1, 2.1, 3.1], name = "testnm", group = "testgrp"},
                {m = 2.0, pos = [1.2, 2.2, 3.2], vel = [1.3, 2.3, 3.3], name = "testnm", group = "testgrp"},
                {m = 3.0, pos = [1.4, 2.4, 3.4], vel = [1.5, 2.5, 3.5]}
            ]
        )"_toml;

        const auto sys = mjolnir::read_system<traits_type>(v, 0);

        BOOST_TEST(sys.boundary().lower_bound().at(0) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.boundary().lower_bound().at(1) == 1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.boundary().lower_bound().at(2) == 2.0, boost::test_tools::tolerance(tol));

        BOOST_TEST(sys.boundary().upper_bound().at(0) == 10.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.boundary().upper_bound().at(1) == 11.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.boundary().upper_bound().at(2) == 12.0, boost::test_tools::tolerance(tol));

        BOOST_TEST(sys.size() == 3u);

        BOOST_TEST(sys.mass(0) == 1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.mass(1) == 2.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.mass(2) == 3.0, boost::test_tools::tolerance(tol));

        BOOST_TEST(sys.rmass(0) == 1.0 / 1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.rmass(1) == 1.0 / 2.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.rmass(2) == 1.0 / 3.0, boost::test_tools::tolerance(tol));

        BOOST_TEST(sys.position(0).at(0) == 1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.position(0).at(1) == 2.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.position(0).at(2) == 3.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.position(1).at(0) == 1.2, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.position(1).at(1) == 2.2, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.position(1).at(2) == 3.2, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.position(2).at(0) == 1.4, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.position(2).at(1) == 2.4, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.position(2).at(2) == 3.4, boost::test_tools::tolerance(tol));

        BOOST_TEST(sys.velocity(0).at(0) == 1.1, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.velocity(0).at(1) == 2.1, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.velocity(0).at(2) == 3.1, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.velocity(1).at(0) == 1.3, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.velocity(1).at(1) == 2.3, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.velocity(1).at(2) == 3.3, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.velocity(2).at(0) == 1.5, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.velocity(2).at(1) == 2.5, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.velocity(2).at(2) == 3.5, boost::test_tools::tolerance(tol));

        BOOST_TEST(sys.force(0).at(0) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(0).at(1) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(0).at(2) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(1).at(0) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(1).at(1) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(1).at(2) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(2).at(0) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(2).at(1) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys.force(2).at(2) == 0.0, boost::test_tools::tolerance(tol));

        BOOST_TEST(sys.name(0) == "testnm");
        BOOST_TEST(sys.name(1) == "testnm");
        BOOST_TEST(sys.name(2) == "X");

        BOOST_TEST(sys.group(0) == "testgrp");
        BOOST_TEST(sys.group(1) == "testgrp");
        BOOST_TEST(sys.group(2) == "NONE");

        BOOST_TEST(sys.has_attribute("test_attr"));
        BOOST_TEST(sys.attribute("test_attr") == 3.14, boost::test_tools::tolerance(tol));
    }
}
