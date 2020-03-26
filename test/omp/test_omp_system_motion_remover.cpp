#define BOOST_TEST_MODULE "test_omp_system_motion_remover"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SystemMotionRemover.hpp>
#include <mjolnir/omp/SystemMotionRemover.hpp>
#include <random>

BOOST_AUTO_TEST_CASE(SystemMotionRemover_both)
{
    mjolnir::LoggerManager::set_default_logger("test_omp_system_motion_remover.log");

    using traits_type   = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type     = traits_type::real_type;
    using coord_type    = traits_type::coordinate_type;
    using boundary_type = traits_type::boundary_type;
    using system_type   = mjolnir::System<traits_type>;

    constexpr real_type tol = 1e-8;
    constexpr std::size_t N = 100;

    system_type sys(N, boundary_type{});
    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(0.0, 1.0);

    for(std::size_t n=0; n<100; ++n)
    {
        real_type E_kinetic_pre = 0.0;
        auto CoM = mjolnir::math::make_coordinate<coord_type>(0,0,0);
        for(std::size_t i=0; i<N; ++i)
        {
            sys.mass(i)     = 1.0 + uni(mt);
            sys.rmass(i)    = 1.0 / sys.mass(i);
            sys.position(i) = mjolnir::math::make_coordinate<coord_type>(uni(mt), uni(mt), uni(mt));
            sys.velocity(i) = mjolnir::math::make_coordinate<coord_type>(uni(mt), uni(mt), uni(mt));
            sys.force(i)    = mjolnir::math::make_coordinate<coord_type>(0,0,0);

            sys.name(0)  = "X";
            sys.group(0) = "X";

            CoM += sys.mass(i) * sys.position(i);
            E_kinetic_pre += 0.5 * sys.mass(i) * mjolnir::math::length_sq(sys.velocity(i));
        }
        const system_type sys_pre(sys);

        mjolnir::SystemMotionRemover<traits_type> remover(true, true, true);

        remover.remove(sys);

        auto translation = mjolnir::math::make_coordinate<coord_type>(0,0,0);
        auto rotation    = mjolnir::math::make_coordinate<coord_type>(0,0,0);

        real_type E_kinetic_post = 0.0;
        for(std::size_t i=0; i<N; ++i)
        {
            const auto  m = sys.mass(i);
            const auto  r = sys.position(i) - CoM;
            const auto& v = sys.velocity(i);

            translation += m * v;
            rotation    += m * mjolnir::math::cross_product(r, v);

            E_kinetic_post += 0.5 * sys.mass(i) * mjolnir::math::length_sq(sys.velocity(i));

            BOOST_TEST(mjolnir::math::X(sys.position(i)) == mjolnir::math::X(sys_pre.position(i)));
            BOOST_TEST(mjolnir::math::Y(sys.position(i)) == mjolnir::math::Y(sys_pre.position(i)));
            BOOST_TEST(mjolnir::math::Z(sys.position(i)) == mjolnir::math::Z(sys_pre.position(i)));
        }

        BOOST_TEST(mjolnir::math::X(translation) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Y(translation) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Z(translation) == 0.0, boost::test_tools::tolerance(tol));

        BOOST_TEST(mjolnir::math::X(rotation) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Y(rotation) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Z(rotation) == 0.0, boost::test_tools::tolerance(tol));

        // rescaled.
        BOOST_TEST(E_kinetic_pre == E_kinetic_post, boost::test_tools::tolerance(tol));
    }

    for(std::size_t n=0; n<100; ++n)
    {
        real_type E_kinetic_pre = 0.0;
        auto CoM = mjolnir::math::make_coordinate<coord_type>(0,0,0);
        for(std::size_t i=0; i<N; ++i)
        {
            sys.mass(i)     = 1.0 + uni(mt);
            sys.rmass(i)    = 1.0 / sys.mass(i);
            sys.position(i) = mjolnir::math::make_coordinate<coord_type>(uni(mt), uni(mt), uni(mt));
            sys.velocity(i) = mjolnir::math::make_coordinate<coord_type>(uni(mt), uni(mt), uni(mt));
            sys.force(i)    = mjolnir::math::make_coordinate<coord_type>(0,0,0);

            sys.name(0)  = "X";
            sys.group(0) = "X";

            CoM += sys.mass(i) * sys.position(i);
            E_kinetic_pre += 0.5 * sys.mass(i) * mjolnir::math::length_sq(sys.velocity(i));
        }

        const system_type sys_pre(sys);

        mjolnir::SystemMotionRemover<traits_type> remover(true, true, false);

        remover.remove(sys);

        auto translation = mjolnir::math::make_coordinate<coord_type>(0,0,0);
        auto rotation    = mjolnir::math::make_coordinate<coord_type>(0,0,0);

        real_type E_kinetic_post = 0.0;
        for(std::size_t i=0; i<N; ++i)
        {
            const auto  m = sys.mass(i);
            const auto  r = sys.position(i) - CoM;
            const auto& v = sys.velocity(i);

            translation += m * v;
            rotation    += m * mjolnir::math::cross_product(r, v);

            E_kinetic_post += 0.5 * sys.mass(i) * mjolnir::math::length_sq(sys.velocity(i));

            BOOST_TEST(mjolnir::math::X(sys.position(i)) == mjolnir::math::X(sys_pre.position(i)));
            BOOST_TEST(mjolnir::math::Y(sys.position(i)) == mjolnir::math::Y(sys_pre.position(i)));
            BOOST_TEST(mjolnir::math::Z(sys.position(i)) == mjolnir::math::Z(sys_pre.position(i)));
        }

        BOOST_TEST(mjolnir::math::X(translation) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Y(translation) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Z(translation) == 0.0, boost::test_tools::tolerance(tol));

        BOOST_TEST(mjolnir::math::X(rotation) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Y(rotation) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Z(rotation) == 0.0, boost::test_tools::tolerance(tol));

        // not rescaled.
        BOOST_TEST(E_kinetic_pre >= E_kinetic_post);
    }
}

BOOST_AUTO_TEST_CASE(SystemMotionRemover_translation)
{
    mjolnir::LoggerManager::set_default_logger("test_omp_system_motion_remover.log");

    using traits_type   = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type     = traits_type::real_type;
    using coord_type    = traits_type::coordinate_type;
    using boundary_type = traits_type::boundary_type;
    using system_type   = mjolnir::System<traits_type>;

    constexpr real_type tol = 1e-8;
    constexpr std::size_t N = 100;

    system_type sys(N, boundary_type{});
    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(0.0, 1.0);

    for(std::size_t n=0; n<100; ++n)
    {
        real_type E_kinetic_pre = 0.0;
        for(std::size_t i=0; i<N; ++i)
        {
            sys.mass(i)     = 1.0 + uni(mt);
            sys.rmass(i)    = 1.0 / sys.mass(i);
            sys.position(i) = mjolnir::math::make_coordinate<coord_type>(uni(mt), uni(mt), uni(mt));
            sys.velocity(i) = mjolnir::math::make_coordinate<coord_type>(uni(mt), uni(mt), uni(mt));
            sys.force(i)    = mjolnir::math::make_coordinate<coord_type>(0,0,0);

            sys.name(0)  = "X";
            sys.group(0) = "X";

            E_kinetic_pre += 0.5 * sys.mass(i) * mjolnir::math::length_sq(sys.velocity(i));
        }
        const system_type sys_pre(sys);

        mjolnir::SystemMotionRemover<traits_type> remover(true, false, true);

        remover.remove(sys);

        auto translation = mjolnir::math::make_coordinate<coord_type>(0,0,0);

        real_type E_kinetic_post = 0.0;
        for(std::size_t i=0; i<N; ++i)
        {
            const auto  m = sys.mass(i);
            const auto& v = sys.velocity(i);

            translation += m * v;
            E_kinetic_post += 0.5 * m * mjolnir::math::length_sq(v);

            BOOST_TEST(mjolnir::math::X(sys.position(i)) == mjolnir::math::X(sys_pre.position(i)));
            BOOST_TEST(mjolnir::math::Y(sys.position(i)) == mjolnir::math::Y(sys_pre.position(i)));
            BOOST_TEST(mjolnir::math::Z(sys.position(i)) == mjolnir::math::Z(sys_pre.position(i)));
        }

        BOOST_TEST(mjolnir::math::X(translation) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Y(translation) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Z(translation) == 0.0, boost::test_tools::tolerance(tol));

        // rescaled.
        BOOST_TEST(E_kinetic_pre == E_kinetic_post, boost::test_tools::tolerance(tol));
    }

    for(std::size_t n=0; n<100; ++n)
    {
        real_type E_kinetic_pre = 0.0;
        for(std::size_t i=0; i<N; ++i)
        {
            sys.mass(i)     = 1.0 + uni(mt);
            sys.rmass(i)    = 1.0 / sys.mass(i);
            sys.position(i) = mjolnir::math::make_coordinate<coord_type>(uni(mt), uni(mt), uni(mt));
            sys.velocity(i) = mjolnir::math::make_coordinate<coord_type>(uni(mt), uni(mt), uni(mt));
            sys.force(i)    = mjolnir::math::make_coordinate<coord_type>(0,0,0);

            sys.name(0)  = "X";
            sys.group(0) = "X";

            E_kinetic_pre += 0.5 * sys.mass(i) * mjolnir::math::length_sq(sys.velocity(i));
        }
        const system_type sys_pre(sys);
        mjolnir::SystemMotionRemover<traits_type> remover(true, false, false);

        remover.remove(sys);

        auto translation = mjolnir::math::make_coordinate<coord_type>(0,0,0);

        real_type E_kinetic_post = 0.0;
        for(std::size_t i=0; i<N; ++i)
        {
            const auto  m = sys.mass(i);
            const auto& v = sys.velocity(i);

            translation += m * v;

            E_kinetic_post += 0.5 * m * mjolnir::math::length_sq(v);

            BOOST_TEST(mjolnir::math::X(sys.position(i)) == mjolnir::math::X(sys_pre.position(i)));
            BOOST_TEST(mjolnir::math::Y(sys.position(i)) == mjolnir::math::Y(sys_pre.position(i)));
            BOOST_TEST(mjolnir::math::Z(sys.position(i)) == mjolnir::math::Z(sys_pre.position(i)));
        }

        BOOST_TEST(mjolnir::math::X(translation) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Y(translation) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Z(translation) == 0.0, boost::test_tools::tolerance(tol));

        // not rescaled.
        BOOST_TEST(E_kinetic_pre >= E_kinetic_post);
    }
}

BOOST_AUTO_TEST_CASE(SystemMotionRemover_rotation)
{
    mjolnir::LoggerManager::set_default_logger("test_omp_system_motion_remover.log");

    using traits_type   = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type     = traits_type::real_type;
    using coord_type    = traits_type::coordinate_type;
    using boundary_type = traits_type::boundary_type;
    using system_type   = mjolnir::System<traits_type>;

    constexpr real_type tol = 1e-8;
    constexpr std::size_t N = 100;

    system_type sys(N, boundary_type{});
    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(0.0, 1.0);

    for(std::size_t n=0; n<100; ++n)
    {
        real_type E_kinetic_pre = 0.0;
        auto CoM = mjolnir::math::make_coordinate<coord_type>(0,0,0);
        for(std::size_t i=0; i<N; ++i)
        {
            sys.mass(i)     = 1.0 + uni(mt);
            sys.rmass(i)    = 1.0 / sys.mass(i);
            sys.position(i) = mjolnir::math::make_coordinate<coord_type>(uni(mt), uni(mt), uni(mt));
            sys.velocity(i) = mjolnir::math::make_coordinate<coord_type>(uni(mt), uni(mt), uni(mt));
            sys.force(i)    = mjolnir::math::make_coordinate<coord_type>(0,0,0);

            sys.name(0)  = "X";
            sys.group(0) = "X";

            CoM += sys.mass(i) * sys.position(i);
            E_kinetic_pre += 0.5 * sys.mass(i) * mjolnir::math::length_sq(sys.velocity(i));
        }
        const system_type sys_pre(sys);
        mjolnir::SystemMotionRemover<traits_type> remover(true, true, true);

        remover.remove(sys);

        auto rotation    = mjolnir::math::make_coordinate<coord_type>(0,0,0);

        real_type E_kinetic_post = 0.0;
        for(std::size_t i=0; i<N; ++i)
        {
            const auto  m = sys.mass(i);
            const auto  r = sys.position(i) - CoM;
            const auto& v = sys.velocity(i);

            rotation       += m * mjolnir::math::cross_product(r, v);
            E_kinetic_post += 0.5 * m * mjolnir::math::length_sq(v);

            BOOST_TEST(mjolnir::math::X(sys.position(i)) == mjolnir::math::X(sys_pre.position(i)));
            BOOST_TEST(mjolnir::math::Y(sys.position(i)) == mjolnir::math::Y(sys_pre.position(i)));
            BOOST_TEST(mjolnir::math::Z(sys.position(i)) == mjolnir::math::Z(sys_pre.position(i)));
        }

        BOOST_TEST(mjolnir::math::X(rotation) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Y(rotation) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Z(rotation) == 0.0, boost::test_tools::tolerance(tol));

        // rescaled.
        BOOST_TEST(E_kinetic_pre == E_kinetic_post, boost::test_tools::tolerance(tol));
    }

    for(std::size_t n=0; n<100; ++n)
    {
        real_type E_kinetic_pre = 0.0;
        auto CoM = mjolnir::math::make_coordinate<coord_type>(0,0,0);
        for(std::size_t i=0; i<N; ++i)
        {
            sys.mass(i)     = 1.0 + uni(mt);
            sys.rmass(i)    = 1.0 / sys.mass(i);
            sys.position(i) = mjolnir::math::make_coordinate<coord_type>(uni(mt), uni(mt), uni(mt));
            sys.velocity(i) = mjolnir::math::make_coordinate<coord_type>(uni(mt), uni(mt), uni(mt));
            sys.force(i)    = mjolnir::math::make_coordinate<coord_type>(0,0,0);

            sys.name(0)  = "X";
            sys.group(0) = "X";

            CoM += sys.mass(i) * sys.position(i);
            E_kinetic_pre += 0.5 * sys.mass(i) * mjolnir::math::length_sq(sys.velocity(i));
        }
        const system_type sys_pre(sys);
        mjolnir::SystemMotionRemover<traits_type> remover(true, true, false);

        remover.remove(sys);

        auto rotation    = mjolnir::math::make_coordinate<coord_type>(0,0,0);

        real_type E_kinetic_post = 0.0;
        for(std::size_t i=0; i<N; ++i)
        {
            const auto  m = sys.mass(i);
            const auto  r = sys.position(i) - CoM;
            const auto& v = sys.velocity(i);

            rotation       += m * mjolnir::math::cross_product(r, v);
            E_kinetic_post += 0.5 * m * mjolnir::math::length_sq(v);

            BOOST_TEST(mjolnir::math::X(sys.position(i)) == mjolnir::math::X(sys_pre.position(i)));
            BOOST_TEST(mjolnir::math::Y(sys.position(i)) == mjolnir::math::Y(sys_pre.position(i)));
            BOOST_TEST(mjolnir::math::Z(sys.position(i)) == mjolnir::math::Z(sys_pre.position(i)));
        }

        BOOST_TEST(mjolnir::math::X(rotation) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Y(rotation) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Z(rotation) == 0.0, boost::test_tools::tolerance(tol));

        // not rescaled.
        BOOST_TEST(E_kinetic_pre >= E_kinetic_post);
    }
}
