#define BOOST_TEST_MODULE "test_omp_bond_angle_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/math/math.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/MultipleBasin/MultipleBasinForceField.hpp>
#include <mjolnir/forcefield/MultipleBasin/MultipleBasin2BasinUnit.hpp>
#include <mjolnir/forcefield/MultipleBasin/MultipleBasin3BasinUnit.hpp>
#include <mjolnir/forcefield/local/HarmonicPotential.hpp>
#include <mjolnir/omp/BondLengthInteraction.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/RandomNumberGenerator.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(omp_MultipleBasin_2Basin_numerical_difference)
{
    namespace test = mjolnir::test;

    mjolnir::LoggerManager::set_default_logger("test_omp_multiple_basin_forcefield.log");
    using traits_type      = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::HarmonicPotential<real_type>;
    using interaction_type = mjolnir::BondLengthInteraction<traits_type, potential_type>;

    const real_type k(10.0);
    const real_type native1(1.5);
    const real_type native2(2.0);

    mjolnir::LocalForceField<traits_type>    loc1, loc2, locc;
    mjolnir::GlobalForceField<traits_type>   glo1, glo2, gloc;
    mjolnir::ExternalForceField<traits_type> ext1, ext2, extc;

    std::vector<std::pair<std::array<std::size_t, 2>, potential_type>> param1, param2;
    param1.emplace_back(std::array<std::size_t, 2>{{0,1}}, potential_type(k, native1));
    param2.emplace_back(std::array<std::size_t, 2>{{0,1}}, potential_type(k, native2));

    loc1.emplace(mjolnir::make_unique<interaction_type>("none", std::move(param1)));
    loc2.emplace(mjolnir::make_unique<interaction_type>("none", std::move(param2)));

    std::unique_ptr<mjolnir::MultipleBasinUnitBase<traits_type>> unit1 =
        mjolnir::make_unique<mjolnir::MultipleBasin2BasinUnit<traits_type>>(
            -10.0, "short", "long", 0.0, 0.0,
            std::make_tuple(std::move(loc1), std::move(glo1), std::move(ext1)),
            std::make_tuple(std::move(loc2), std::move(glo2), std::move(ext2)),
            0.0, 0.0);

    std::vector<std::unique_ptr<mjolnir::MultipleBasinUnitBase<traits_type>>> units;
    units.push_back(std::move(unit1));

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    mjolnir::MultipleBasinForceField<traits_type> forcefield(
            std::make_tuple(std::move(locc), std::move(gloc), std::move(extc)),
            std::move(units));

    for(std::size_t i = 0; i < 100; ++i)
    {
        system_type sys(2, boundary_type{});

        test::clear_everything(sys);

        sys.position(0) = coord_type( 0.0 + 0.01 * uni(mt), 0.0 + 0.01 * uni(mt), 0.0 + 0.01 * uni(mt));
        sys.position(1) = coord_type( 1.0 + 0.01 * uni(mt), 1.0 + 0.01 * uni(mt), 1.0 + 0.01 * uni(mt));

        forcefield.initialize(sys); // don't forget this

        constexpr real_type tol = 1e-4;
        constexpr real_type dr  = 1e-5;

        test::check_force(sys, forcefield, tol, dr);
    }
}

// check difference of forces
BOOST_AUTO_TEST_CASE(omp_MultipleBasin_2Basin_consistency)
{
    namespace test = mjolnir::test;

    mjolnir::LoggerManager::set_default_logger("test_omp_multiple_basin_forcefield.log");
    using traits_type      = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::HarmonicPotential<real_type>;
    using interaction_type = mjolnir::BondLengthInteraction<traits_type, potential_type>;

    using default_traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using default_system_type      = mjolnir::System<default_traits_type>;
    using default_potential_type   = mjolnir::HarmonicPotential<real_type>;
    using default_interaction_type = mjolnir::BondLengthInteraction<default_traits_type, default_potential_type>;

    const real_type k(10.0);
    const real_type native1(1.5);
    const real_type native2(2.0);
    const std::size_t Nparticle = 1000;

    mjolnir::MultipleBasinForceField<traits_type> ff_openmp = [&]()
    {
        mjolnir::LocalForceField<traits_type>    loc1, loc2, locc;
        mjolnir::GlobalForceField<traits_type>   glo1, glo2, gloc;
        mjolnir::ExternalForceField<traits_type> ext1, ext2, extc;

        std::vector<std::pair<std::array<std::size_t, 2>, potential_type>> param1, param2;
        for(std::size_t i=0; i<Nparticle-1; ++i)
        {
            param1.emplace_back(std::array<std::size_t, 2>{{i,i+1}}, potential_type(k, native1));
            param2.emplace_back(std::array<std::size_t, 2>{{i,i+1}}, potential_type(k, native2));
        }
        loc1.emplace(mjolnir::make_unique<interaction_type>("none", std::move(param1)));
        loc2.emplace(mjolnir::make_unique<interaction_type>("none", std::move(param2)));

        std::unique_ptr<mjolnir::MultipleBasinUnitBase<traits_type>> unit1 =
            mjolnir::make_unique<mjolnir::MultipleBasin2BasinUnit<traits_type>>(
                -10.0, "short", "long", 0.0, 0.0,
                std::make_tuple(std::move(loc1), std::move(glo1), std::move(ext1)),
                std::make_tuple(std::move(loc2), std::move(glo2), std::move(ext2)),
                0.0, 0.0);

        std::vector<std::unique_ptr<mjolnir::MultipleBasinUnitBase<traits_type>>> units;
        units.push_back(std::move(unit1));

        return mjolnir::MultipleBasinForceField<traits_type>(
                std::make_tuple(std::move(locc), std::move(gloc), std::move(extc)),
                std::move(units));
    }();

    mjolnir::MultipleBasinForceField<default_traits_type> ff_default = [&]()
    {
        mjolnir::LocalForceField<default_traits_type>    loc1, loc2, locc;
        mjolnir::GlobalForceField<default_traits_type>   glo1, glo2, gloc;
        mjolnir::ExternalForceField<default_traits_type> ext1, ext2, extc;

        std::vector<std::pair<std::array<std::size_t, 2>, default_potential_type>> param1, param2;
        for(std::size_t i=0; i<Nparticle-1; ++i)
        {
            param1.emplace_back(std::array<std::size_t, 2>{{i,i+1}}, default_potential_type(k, native1));
            param2.emplace_back(std::array<std::size_t, 2>{{i,i+1}}, default_potential_type(k, native2));
        }
        loc1.emplace(mjolnir::make_unique<default_interaction_type>("none", std::move(param1)));
        loc2.emplace(mjolnir::make_unique<default_interaction_type>("none", std::move(param2)));

        std::unique_ptr<mjolnir::MultipleBasinUnitBase<default_traits_type>> unit1 =
            mjolnir::make_unique<mjolnir::MultipleBasin2BasinUnit<default_traits_type>>(
                -10.0, "short", "long", 0.0, 0.0,
                std::make_tuple(std::move(loc1), std::move(glo1), std::move(ext1)),
                std::make_tuple(std::move(loc2), std::move(glo2), std::move(ext2)),
                0.0, 0.0);

        std::vector<std::unique_ptr<mjolnir::MultipleBasinUnitBase<default_traits_type>>> units;
        units.push_back(std::move(unit1));

        return mjolnir::MultipleBasinForceField<default_traits_type>(
                std::make_tuple(std::move(locc), std::move(gloc), std::move(extc)),
                std::move(units));
    }();

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    for(std::size_t i = 0; i < 100; ++i)
    {
        system_type         sys_openmp (Nparticle, boundary_type{});
        default_system_type sys_default(Nparticle, boundary_type{});

        test::clear_everything(sys_openmp);
        test::clear_everything(sys_default);

        for(std::size_t j=0; j<Nparticle; ++j)
        {
            const coord_type pos(j * 1.0 + 0.01 * uni(mt),
                                 j * 1.0 + 0.01 * uni(mt),
                                 j * 1.0 + 0.01 * uni(mt));

            sys_openmp.position(j) = pos;
            sys_default.position(j) = pos;
        }

        ff_openmp .initialize(sys_openmp ); // don't forget this
        ff_default.initialize(sys_default); // don't forget this

        constexpr real_type tol = 1e-4;

        // --------------------------------------------------------------------
        // check_force

        sys_openmp.preprocess_forces();
        ff_openmp .calc_force(sys_openmp);
        sys_openmp.postprocess_forces();

        ff_default.calc_force(sys_default);

        for(std::size_t j=0; j<Nparticle; ++j)
        {
            using namespace mjolnir::math;
            BOOST_TEST(X(sys_default.force(j)) == X(sys_openmp.force(j)),
                       boost::test_tools::tolerance(tol));

            BOOST_TEST(Y(sys_default.force(j)) == Y(sys_openmp.force(j)),
                       boost::test_tools::tolerance(tol));

            BOOST_TEST(Z(sys_default.force(j)) == Z(sys_openmp.force(j)),
                       boost::test_tools::tolerance(tol));
        }
        // check the virials are the same
        for(std::size_t j=0; j<9; ++j)
        {
            BOOST_TEST(sys_default.virial()[j] == sys_openmp.virial()[j], boost::test_tools::tolerance(tol));
        }

        // --------------------------------------------------------------------
        // check_force_and_virial

        test::clear_force(sys_openmp);
        test::clear_force(sys_default);

        sys_openmp.preprocess_forces();
        ff_openmp .calc_force_and_virial(sys_openmp);
        sys_openmp.postprocess_forces();

        ff_default.calc_force_and_virial(sys_default);

        for(std::size_t j=0; j<Nparticle; ++j)
        {
            using namespace mjolnir::math;
            BOOST_TEST(X(sys_default.force(j)) == X(sys_openmp.force(j)),
                       boost::test_tools::tolerance(tol));

            BOOST_TEST(Y(sys_default.force(j)) == Y(sys_openmp.force(j)),
                       boost::test_tools::tolerance(tol));

            BOOST_TEST(Z(sys_default.force(j)) == Z(sys_openmp.force(j)),
                       boost::test_tools::tolerance(tol));
        }
        // check the virials are the same
        for(std::size_t j=0; j<9; ++j)
        {
            BOOST_TEST(sys_default.virial()[j] == sys_openmp.virial()[j], boost::test_tools::tolerance(tol));
        }

        // --------------------------------------------------------------------
        // check_energy

        test::clear_force(sys_openmp);
        test::clear_force(sys_default);

        const auto omp_ene = ff_openmp .calc_energy(sys_openmp);
        const auto seq_ene = ff_default.calc_energy(sys_default);

        BOOST_TEST(omp_ene == seq_ene, boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(MultipleBasin_3Basin_numerical_difference)
{
    namespace test = mjolnir::test;

    mjolnir::LoggerManager::set_default_logger("test_omp_multiple_basin_forcefield.log");
    using traits_type      = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::HarmonicPotential<real_type>;
    using interaction_type = mjolnir::BondLengthInteraction<traits_type, potential_type>;

    const real_type k(10.0);
    const real_type native1(1.6);
    const real_type native2(1.8);
    const real_type native3(1.9);

    mjolnir::LocalForceField<traits_type>    loc1, loc2, loc3, locc;
    mjolnir::GlobalForceField<traits_type>   glo1, glo2, glo3, gloc;
    mjolnir::ExternalForceField<traits_type> ext1, ext2, ext3, extc;

    std::vector<std::pair<std::array<std::size_t, 2>, potential_type>> param1, param2, param3;
    param1.emplace_back(std::array<std::size_t, 2>{{0,1}}, potential_type(k, native1));
    param2.emplace_back(std::array<std::size_t, 2>{{0,1}}, potential_type(k, native2));
    param3.emplace_back(std::array<std::size_t, 2>{{0,1}}, potential_type(k, native3));

    loc1.emplace(mjolnir::make_unique<interaction_type>("none", std::move(param1)));
    loc2.emplace(mjolnir::make_unique<interaction_type>("none", std::move(param2)));
    loc3.emplace(mjolnir::make_unique<interaction_type>("none", std::move(param3)));

    std::unique_ptr<mjolnir::MultipleBasinUnitBase<traits_type>> unit1 =
        mjolnir::make_unique<mjolnir::MultipleBasin3BasinUnit<traits_type>>(
            "short", "middle", "long", -10.0, -10.0, -10.0, 0.0, 0.0, 0.0,
            std::make_tuple(std::move(loc1), std::move(glo1), std::move(ext1)),
            std::make_tuple(std::move(loc2), std::move(glo2), std::move(ext2)),
            std::make_tuple(std::move(loc3), std::move(glo3), std::move(ext3)));

    std::vector<std::unique_ptr<mjolnir::MultipleBasinUnitBase<traits_type>>> units;
    units.push_back(std::move(unit1));

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    mjolnir::MultipleBasinForceField<traits_type> forcefield(
            std::make_tuple(std::move(locc), std::move(gloc), std::move(extc)),
            std::move(units));

    for(std::size_t i = 0; i < 100; ++i)
    {
        system_type sys(2, boundary_type{});

        test::clear_everything(sys);

        sys.position(0) = coord_type( 0.0 + 0.01 * uni(mt), 0.0 + 0.01 * uni(mt), 0.0 + 0.01 * uni(mt));
        sys.position(1) = coord_type( 1.0 + 0.01 * uni(mt), 1.0 + 0.01 * uni(mt), 1.0 + 0.01 * uni(mt));

        forcefield.initialize(sys); // don't forget this

        constexpr real_type tol = 1e-4;
        constexpr real_type dr  = 1e-5;
        test::check_force(sys, forcefield, tol, dr);
    }
}

// check difference of forces
BOOST_AUTO_TEST_CASE(omp_MultipleBasin_3Basin_consistency)
{
    namespace test = mjolnir::test;

    mjolnir::LoggerManager::set_default_logger("test_omp_multiple_basin_forcefield.log");
    using traits_type      = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::HarmonicPotential<real_type>;
    using interaction_type = mjolnir::BondLengthInteraction<traits_type, potential_type>;

    using default_traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using default_system_type      = mjolnir::System<default_traits_type>;
    using default_potential_type   = mjolnir::HarmonicPotential<real_type>;
    using default_interaction_type = mjolnir::BondLengthInteraction<default_traits_type, default_potential_type>;

    const real_type k(10.0);
    const real_type native1(1.6);
    const real_type native2(1.8);
    const real_type native3(1.9);

    const std::size_t Nparticle = 1000;

    mjolnir::MultipleBasinForceField<traits_type> ff_openmp = [&]()
    {
        mjolnir::LocalForceField<traits_type>    loc1, loc2, loc3, locc;
        mjolnir::GlobalForceField<traits_type>   glo1, glo2, glo3, gloc;
        mjolnir::ExternalForceField<traits_type> ext1, ext2, ext3, extc;

        std::vector<std::pair<std::array<std::size_t, 2>, potential_type>> param1, param2, param3;
        for(std::size_t i=0; i<Nparticle-1; ++i)
        {
            param1.emplace_back(std::array<std::size_t, 2>{{i,i+1}}, potential_type(k, native1));
            param2.emplace_back(std::array<std::size_t, 2>{{i,i+1}}, potential_type(k, native2));
            param3.emplace_back(std::array<std::size_t, 2>{{i,i+1}}, potential_type(k, native3));
        }
        loc1.emplace(mjolnir::make_unique<interaction_type>("none", std::move(param1)));
        loc2.emplace(mjolnir::make_unique<interaction_type>("none", std::move(param2)));
        loc3.emplace(mjolnir::make_unique<interaction_type>("none", std::move(param3)));

        std::unique_ptr<mjolnir::MultipleBasinUnitBase<traits_type>> unit1 =
            mjolnir::make_unique<mjolnir::MultipleBasin3BasinUnit<traits_type>>(
                "short", "middle", "long", -10.0, -10.0, -10.0, 0.0, 0.0, 0.0,
                std::make_tuple(std::move(loc1), std::move(glo1), std::move(ext1)),
                std::make_tuple(std::move(loc2), std::move(glo2), std::move(ext2)),
                std::make_tuple(std::move(loc3), std::move(glo3), std::move(ext3)));

        std::vector<std::unique_ptr<mjolnir::MultipleBasinUnitBase<traits_type>>> units;
        units.push_back(std::move(unit1));

        return mjolnir::MultipleBasinForceField<traits_type>(
                std::make_tuple(std::move(locc), std::move(gloc), std::move(extc)),
                std::move(units));
    }();

    mjolnir::MultipleBasinForceField<default_traits_type> ff_default = [&]()
    {
        mjolnir::LocalForceField<default_traits_type>    loc1, loc2, loc3, locc;
        mjolnir::GlobalForceField<default_traits_type>   glo1, glo2, glo3, gloc;
        mjolnir::ExternalForceField<default_traits_type> ext1, ext2, ext3, extc;

        std::vector<std::pair<std::array<std::size_t, 2>, default_potential_type>> param1, param2, param3;
        for(std::size_t i=0; i<Nparticle-1; ++i)
        {
            param1.emplace_back(std::array<std::size_t, 2>{{i,i+1}}, default_potential_type(k, native1));
            param2.emplace_back(std::array<std::size_t, 2>{{i,i+1}}, default_potential_type(k, native2));
            param3.emplace_back(std::array<std::size_t, 2>{{i,i+1}}, default_potential_type(k, native3));
        }
        loc1.emplace(mjolnir::make_unique<default_interaction_type>("none", std::move(param1)));
        loc2.emplace(mjolnir::make_unique<default_interaction_type>("none", std::move(param2)));
        loc3.emplace(mjolnir::make_unique<default_interaction_type>("none", std::move(param3)));

        std::unique_ptr<mjolnir::MultipleBasinUnitBase<default_traits_type>> unit1 =
            mjolnir::make_unique<mjolnir::MultipleBasin3BasinUnit<default_traits_type>>(
                "short", "middle", "long", -10.0, -10.0, -10.0, 0.0, 0.0, 0.0,
                std::make_tuple(std::move(loc1), std::move(glo1), std::move(ext1)),
                std::make_tuple(std::move(loc2), std::move(glo2), std::move(ext2)),
                std::make_tuple(std::move(loc3), std::move(glo3), std::move(ext3)));

        std::vector<std::unique_ptr<mjolnir::MultipleBasinUnitBase<default_traits_type>>> units;
        units.push_back(std::move(unit1));

        return mjolnir::MultipleBasinForceField<default_traits_type>(
                std::make_tuple(std::move(locc), std::move(gloc), std::move(extc)),
                std::move(units));
    }();

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    for(std::size_t i = 0; i < 100; ++i)
    {
        system_type         sys_openmp (Nparticle, boundary_type{});
        default_system_type sys_default(Nparticle, boundary_type{});

        test::clear_everything(sys_openmp);
        test::clear_everything(sys_default);

        for(std::size_t j=0; j<Nparticle; ++j)
        {
            const coord_type pos(j * 1.0 + 0.01 * uni(mt),
                                 j * 1.0 + 0.01 * uni(mt),
                                 j * 1.0 + 0.01 * uni(mt));

            sys_openmp .position(j) = pos;
            sys_default.position(j) = pos;
        }

        ff_openmp .initialize(sys_openmp ); // don't forget this
        ff_default.initialize(sys_default); // don't forget this

        constexpr real_type tol = 1e-4;

        // --------------------------------------------------------------------
        // check_force

        sys_openmp.preprocess_forces();
        ff_openmp .calc_force(sys_openmp);
        sys_openmp.postprocess_forces();

        ff_default.calc_force(sys_default);

        for(std::size_t j=0; j<Nparticle; ++j)
        {
            using namespace mjolnir::math;
            BOOST_TEST(X(sys_default.force(j)) == X(sys_openmp.force(j)),
                       boost::test_tools::tolerance(tol));

            BOOST_TEST(Y(sys_default.force(j)) == Y(sys_openmp.force(j)),
                       boost::test_tools::tolerance(tol));

            BOOST_TEST(Z(sys_default.force(j)) == Z(sys_openmp.force(j)),
                       boost::test_tools::tolerance(tol));
        }
        // check the virials are the same
        for(std::size_t j=0; j<9; ++j)
        {
            BOOST_TEST(sys_default.virial()[j] == sys_openmp.virial()[j], boost::test_tools::tolerance(tol));
        }

        // --------------------------------------------------------------------
        // check_force_and_virial

        test::clear_force(sys_openmp);
        test::clear_force(sys_default);

        sys_openmp.preprocess_forces();
        ff_openmp .calc_force_and_virial(sys_openmp);
        sys_openmp.postprocess_forces();

        ff_default.calc_force_and_virial(sys_default);

        for(std::size_t j=0; j<Nparticle; ++j)
        {
            using namespace mjolnir::math;
            BOOST_TEST(X(sys_default.force(j)) == X(sys_openmp.force(j)),
                       boost::test_tools::tolerance(tol));

            BOOST_TEST(Y(sys_default.force(j)) == Y(sys_openmp.force(j)),
                       boost::test_tools::tolerance(tol));

            BOOST_TEST(Z(sys_default.force(j)) == Z(sys_openmp.force(j)),
                       boost::test_tools::tolerance(tol));
        }
        // check the virials are the same
        for(std::size_t j=0; j<9; ++j)
        {
            BOOST_TEST(sys_default.virial()[j] == sys_openmp.virial()[j], boost::test_tools::tolerance(tol));
        }

        // --------------------------------------------------------------------
        // check_energy

        test::clear_force(sys_openmp);
        test::clear_force(sys_default);

        const auto omp_ene = ff_openmp .calc_energy(sys_openmp);
        const auto seq_ene = ff_default.calc_energy(sys_default);

        BOOST_TEST(omp_ene == seq_ene, boost::test_tools::tolerance(tol));
    }
}
