#define BOOST_TEST_MODULE "test_save_load_msgpack"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/core/MsgPackSaver.hpp>
#include <mjolnir/input/read_system.hpp>

using real_types = std::tuple<double, float>;

BOOST_AUTO_TEST_CASE_TEMPLATE(test_save_load_msgpack_unlimited, realT, real_types)
{
    using namespace mjolnir;
    LoggerManager::set_default_logger("test_save_load_msgpack.log");
    using traits_type     = SimulatorTraits<realT, UnlimitedBoundary>;
    using coordinate_type = typename traits_type::coordinate_type;
    using system_type     = System<traits_type>;
    using forcefield_type = std::unique_ptr<ForceFieldBase<traits_type>>;
    using msgsaver_type   = MsgPackSaver<traits_type>;
    using boundary_type   = typename system_type::boundary_type;
    using rng_type = RandomNumberGenerator<traits_type>;

    rng_type rng(123456789);
    system_type sys(100, boundary_type());
    forcefield_type ff;

    sys.attribute("pi") = 3.14;
    sys.attribute("e")  = 2.71;

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.mass(i)     = rng.uniform_real01();
        sys.rmass(i)    = realT(1.0) / sys.mass(i);
        sys.position(i) = math::make_coordinate<coordinate_type>(rng.uniform_real01(), rng.uniform_real01(), rng.uniform_real01());
        sys.velocity(i) = math::make_coordinate<coordinate_type>(rng.uniform_real01(), rng.uniform_real01(), rng.uniform_real01());
        sys.force(i)    = math::make_coordinate<coordinate_type>(rng.uniform_real01(), rng.uniform_real01(), rng.uniform_real01());
        sys.name(i)     = std::string("particle") + std::to_string(i);
        sys.group(i)    = std::string("group")    + std::to_string(i);
    }

    {
        msgsaver_type saver("test_save_load_msgpack_unlimited");
        saver.save(sys);
        saver.save(rng);
    }

    MsgPackLoader<traits_type> loader;
    const rng_type loaded_rng = loader.load_rng("test_save_load_msgpack_unlimited_rng.msg");
    const bool same_rng = (rng == loaded_rng);
    BOOST_TEST(same_rng);

    const system_type loaded_sys = load_system_from_msgpack<traits_type>("test_save_load_msgpack_unlimited_system.msg");

    // It should be bitwise-same.
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        BOOST_TEST(sys.mass(i)     == loaded_sys.mass(i)    );
        BOOST_TEST(sys.rmass(i)    == loaded_sys.rmass(i)   );

        BOOST_TEST(math::X(sys.position(i)) == math::X(loaded_sys.position(i)));
        BOOST_TEST(math::Y(sys.position(i)) == math::Y(loaded_sys.position(i)));
        BOOST_TEST(math::Z(sys.position(i)) == math::Z(loaded_sys.position(i)));

        BOOST_TEST(math::X(sys.velocity(i)) == math::X(loaded_sys.velocity(i)));
        BOOST_TEST(math::Y(sys.velocity(i)) == math::Y(loaded_sys.velocity(i)));
        BOOST_TEST(math::Z(sys.velocity(i)) == math::Z(loaded_sys.velocity(i)));

        BOOST_TEST(math::X(sys.force(i))    == math::X(loaded_sys.force(i))   );
        BOOST_TEST(math::Y(sys.force(i))    == math::Y(loaded_sys.force(i))   );
        BOOST_TEST(math::Z(sys.force(i))    == math::Z(loaded_sys.force(i))   );

        BOOST_TEST(sys.name(i)     == loaded_sys.name(i)    );
        BOOST_TEST(sys.group(i)    == loaded_sys.group(i)   );
    }
    BOOST_TEST(sys.attribute("pi") == loaded_sys.attribute("pi"));
    BOOST_TEST(sys.attribute("e")  == loaded_sys.attribute("e") );


}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_save_load_msgpack_periodic, realT, real_types)
{
    using namespace mjolnir;
    LoggerManager::set_default_logger("test_save_load_msgpack.log");
    using traits_type     = SimulatorTraits<realT, CuboidalPeriodicBoundary>;
    using coordinate_type = typename traits_type::coordinate_type;
    using system_type     = System<traits_type>;
    using forcefield_type = std::unique_ptr<ForceFieldBase<traits_type>>;
    using msgsaver_type   = MsgPackSaver<traits_type>;
    using boundary_type   = typename system_type::boundary_type;
    using rng_type        = RandomNumberGenerator<traits_type>;

    rng_type rng(123456789);
    system_type sys(100, boundary_type(
        math::make_coordinate<coordinate_type>(-1.0, -2.0, -3.0),
        math::make_coordinate<coordinate_type>( 1.0,  2.0,  3.0)));

    forcefield_type ff;

    sys.attribute("pi") = 3.14;
    sys.attribute("e")  = 2.71;

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.mass(i)     = rng.uniform_real01();
        sys.rmass(i)    = realT(1.0) / sys.mass(i);
        sys.position(i) = math::make_coordinate<coordinate_type>(rng.uniform_real01(), rng.uniform_real01(), rng.uniform_real01());
        sys.velocity(i) = math::make_coordinate<coordinate_type>(rng.uniform_real01(), rng.uniform_real01(), rng.uniform_real01());
        sys.force(i)    = math::make_coordinate<coordinate_type>(rng.uniform_real01(), rng.uniform_real01(), rng.uniform_real01());
        sys.name(i)     = std::string("particle") + std::to_string(i);
        sys.group(i)    = std::string("group")    + std::to_string(i);
    }

    {
        msgsaver_type saver("test_save_load_msgpack_periodic");
        saver.save(sys);
        saver.save(rng);
    }

    MsgPackLoader<traits_type> loader;
    const rng_type loaded_rng = loader.load_rng("test_save_load_msgpack_periodic_rng.msg");
    const bool same_rng = (rng == loaded_rng);
    BOOST_TEST(same_rng);

    const system_type loaded_sys = load_system_from_msgpack<traits_type>("test_save_load_msgpack_periodic_system.msg");

    // It should be bitwise-same.

    BOOST_TEST(math::X(sys.boundary().lower_bound()) == math::X(loaded_sys.boundary().lower_bound()));
    BOOST_TEST(math::Y(sys.boundary().lower_bound()) == math::Y(loaded_sys.boundary().lower_bound()));
    BOOST_TEST(math::Z(sys.boundary().lower_bound()) == math::Z(loaded_sys.boundary().lower_bound()));

    BOOST_TEST(math::X(sys.boundary().upper_bound()) == math::X(loaded_sys.boundary().upper_bound()));
    BOOST_TEST(math::Y(sys.boundary().upper_bound()) == math::Y(loaded_sys.boundary().upper_bound()));
    BOOST_TEST(math::Z(sys.boundary().upper_bound()) == math::Z(loaded_sys.boundary().upper_bound()));

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        BOOST_TEST(sys.mass(i)     == loaded_sys.mass(i)    );
        BOOST_TEST(sys.rmass(i)    == loaded_sys.rmass(i)   );

        BOOST_TEST(math::X(sys.position(i)) == math::X(loaded_sys.position(i)));
        BOOST_TEST(math::Y(sys.position(i)) == math::Y(loaded_sys.position(i)));
        BOOST_TEST(math::Z(sys.position(i)) == math::Z(loaded_sys.position(i)));

        BOOST_TEST(math::X(sys.velocity(i)) == math::X(loaded_sys.velocity(i)));
        BOOST_TEST(math::Y(sys.velocity(i)) == math::Y(loaded_sys.velocity(i)));
        BOOST_TEST(math::Z(sys.velocity(i)) == math::Z(loaded_sys.velocity(i)));

        BOOST_TEST(math::X(sys.force(i))    == math::X(loaded_sys.force(i))   );
        BOOST_TEST(math::Y(sys.force(i))    == math::Y(loaded_sys.force(i))   );
        BOOST_TEST(math::Z(sys.force(i))    == math::Z(loaded_sys.force(i))   );

        BOOST_TEST(sys.name(i)     == loaded_sys.name(i)    );
        BOOST_TEST(sys.group(i)    == loaded_sys.group(i)   );
    }
    BOOST_TEST(sys.attribute("pi") == loaded_sys.attribute("pi"));
    BOOST_TEST(sys.attribute("e")  == loaded_sys.attribute("e") );
}
