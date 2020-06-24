#define BOOST_TEST_MODULE "test_save_load_msgpack"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/omp/omp.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/core/MsgPackObserver.hpp>
#include <mjolnir/input/read_system.hpp>

using real_types = std::tuple<double, float>;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_save_load_msgpack_unlimited, realT, real_types)
{
    using namespace mjolnir;
    LoggerManager::set_default_logger("test_omp_save_load_msgpack.log");
    using traits_type     = OpenMPSimulatorTraits<realT, UnlimitedBoundary>;
    using coordinate_type = typename traits_type::coordinate_type;
    using system_type     = System<traits_type>;
    using forcefield_type = std::unique_ptr<ForceFieldBase<traits_type>>;
    using msgpackobs_type = MsgPackObserver<traits_type>;
    using boundary_type   = typename system_type::boundary_type;

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<realT> uni(-10.0, 10.0);
    system_type sys(100, boundary_type());
    forcefield_type ff;

    sys.attribute("pi") = 3.14;
    sys.attribute("e")  = 2.71;

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.mass(i)     = uni(mt);
        sys.rmass(i)    = realT(1.0) / sys.mass(i);
        sys.position(i) = math::make_coordinate<coordinate_type>(uni(mt), uni(mt), uni(mt));
        sys.velocity(i) = math::make_coordinate<coordinate_type>(uni(mt), uni(mt), uni(mt));
        sys.force(i)    = math::make_coordinate<coordinate_type>(uni(mt), uni(mt), uni(mt));
        sys.name(i)     = std::string("particle") + std::to_string(i);
        sys.group(i)    = std::string("group")    + std::to_string(i);
    }

    {
        msgpackobs_type obs("test_omp_save_load_msgpack_unlimited");
        obs.initialize(0, 0.1, sys, ff);
        obs.output(0, 0.1, sys, ff);
    }
    const system_type loaded = load_system_from_msgpack<traits_type>("test_omp_save_load_msgpack_unlimited.msg");

    // It should be bitwise-same.
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        BOOST_TEST(sys.mass(i)     == loaded.mass(i)    );
        BOOST_TEST(sys.rmass(i)    == loaded.rmass(i)   );

        BOOST_TEST(math::X(sys.position(i)) == math::X(loaded.position(i)));
        BOOST_TEST(math::Y(sys.position(i)) == math::Y(loaded.position(i)));
        BOOST_TEST(math::Z(sys.position(i)) == math::Z(loaded.position(i)));

        BOOST_TEST(math::X(sys.velocity(i)) == math::X(loaded.velocity(i)));
        BOOST_TEST(math::Y(sys.velocity(i)) == math::Y(loaded.velocity(i)));
        BOOST_TEST(math::Z(sys.velocity(i)) == math::Z(loaded.velocity(i)));

        BOOST_TEST(math::X(sys.force(i))    == math::X(loaded.force(i))   );
        BOOST_TEST(math::Y(sys.force(i))    == math::Y(loaded.force(i))   );
        BOOST_TEST(math::Z(sys.force(i))    == math::Z(loaded.force(i))   );

        BOOST_TEST(sys.name(i)     == loaded.name(i)    );
        BOOST_TEST(sys.group(i)    == loaded.group(i)   );
    }
    BOOST_TEST(sys.attribute("pi") == loaded.attribute("pi"));
    BOOST_TEST(sys.attribute("e")  == loaded.attribute("e") );
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_save_load_msgpack_periodic, realT, real_types)
{
    using namespace mjolnir;
    LoggerManager::set_default_logger("test_omp_save_load_msgpack.log");
    using traits_type     = OpenMPSimulatorTraits<realT, CuboidalPeriodicBoundary>;
    using coordinate_type = typename traits_type::coordinate_type;
    using system_type     = System<traits_type>;
    using forcefield_type = std::unique_ptr<ForceFieldBase<traits_type>>;
    using msgpackobs_type = MsgPackObserver<traits_type>;
    using boundary_type   = typename system_type::boundary_type;

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<realT> uni(-10.0, 10.0);
    system_type sys(100, boundary_type(
        math::make_coordinate<coordinate_type>(-1.0, -2.0, -3.0),
        math::make_coordinate<coordinate_type>( 1.0,  2.0,  3.0)));

    forcefield_type ff;

    sys.attribute("pi") = 3.14;
    sys.attribute("e")  = 2.71;

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.mass(i)     = uni(mt);
        sys.rmass(i)    = realT(1.0) / sys.mass(i);
        sys.position(i) = math::make_coordinate<coordinate_type>(uni(mt), uni(mt), uni(mt));
        sys.velocity(i) = math::make_coordinate<coordinate_type>(uni(mt), uni(mt), uni(mt));
        sys.force(i)    = math::make_coordinate<coordinate_type>(uni(mt), uni(mt), uni(mt));
        sys.name(i)     = std::string("particle") + std::to_string(i);
        sys.group(i)    = std::string("group")    + std::to_string(i);
    }

    {
        msgpackobs_type obs("test_omp_save_load_msgpack_periodic");
        obs.initialize(0, 0.1, sys, ff);
        obs.output(0, 0.1, sys, ff);
    }
    const system_type loaded = load_system_from_msgpack<traits_type>("test_omp_save_load_msgpack_periodic.msg");

    // It should be bitwise-same.

    BOOST_TEST(math::X(sys.boundary().lower_bound()) == math::X(loaded.boundary().lower_bound()));
    BOOST_TEST(math::Y(sys.boundary().lower_bound()) == math::Y(loaded.boundary().lower_bound()));
    BOOST_TEST(math::Z(sys.boundary().lower_bound()) == math::Z(loaded.boundary().lower_bound()));

    BOOST_TEST(math::X(sys.boundary().upper_bound()) == math::X(loaded.boundary().upper_bound()));
    BOOST_TEST(math::Y(sys.boundary().upper_bound()) == math::Y(loaded.boundary().upper_bound()));
    BOOST_TEST(math::Z(sys.boundary().upper_bound()) == math::Z(loaded.boundary().upper_bound()));

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        BOOST_TEST(sys.mass(i)     == loaded.mass(i)    );
        BOOST_TEST(sys.rmass(i)    == loaded.rmass(i)   );

        BOOST_TEST(math::X(sys.position(i)) == math::X(loaded.position(i)));
        BOOST_TEST(math::Y(sys.position(i)) == math::Y(loaded.position(i)));
        BOOST_TEST(math::Z(sys.position(i)) == math::Z(loaded.position(i)));

        BOOST_TEST(math::X(sys.velocity(i)) == math::X(loaded.velocity(i)));
        BOOST_TEST(math::Y(sys.velocity(i)) == math::Y(loaded.velocity(i)));
        BOOST_TEST(math::Z(sys.velocity(i)) == math::Z(loaded.velocity(i)));

        BOOST_TEST(math::X(sys.force(i))    == math::X(loaded.force(i))   );
        BOOST_TEST(math::Y(sys.force(i))    == math::Y(loaded.force(i))   );
        BOOST_TEST(math::Z(sys.force(i))    == math::Z(loaded.force(i))   );

        BOOST_TEST(sys.name(i)     == loaded.name(i)    );
        BOOST_TEST(sys.group(i)    == loaded.group(i)   );
    }
    BOOST_TEST(sys.attribute("pi") == loaded.attribute("pi"));
    BOOST_TEST(sys.attribute("e")  == loaded.attribute("e") );
}
