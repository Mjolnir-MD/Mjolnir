#define BOOST_TEST_MODULE "test_bond_length_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/MultipleBasin/MultipleBasinForceField.hpp>
#include <mjolnir/forcefield/MultipleBasin/MultipleBasin2BasinUnit.hpp>
#include <mjolnir/forcefield/MultipleBasin/MultipleBasin3BasinUnit.hpp>
#include <mjolnir/forcefield/local/BondLengthInteraction.hpp>
#include <mjolnir/forcefield/local/HarmonicPotential.hpp>
#include <mjolnir/util/make_unique.hpp>

#include <random>

BOOST_AUTO_TEST_CASE(MultipleBasin_2Basin_numerical_difference)
{
    namespace test = mjolnir::test;

    mjolnir::LoggerManager::set_default_logger("test_multiple_basin_forcefield.log");

    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
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

    mjolnir::MultipleBasinForceField<traits_type> forcefield(
            std::make_tuple(std::move(locc), std::move(gloc), std::move(extc)),
            std::move(units));

    for(std::size_t i = 0; i < 1000; ++i)
    {
        system_type sys(2, boundary_type{});

        test::clear_everything(sys);

        sys.position(0) = coord_type( 0.0, 0.0, 0.0);
        sys.position(1) = coord_type( 1.0, 1.0, 1.0);

        test::apply_random_rotation(sys, mt);
        test::apply_random_perturbation(sys, mt, 0.01);

        forcefield.initialize(sys); // don't forget this

        constexpr real_type tol = 1e-4;
        constexpr real_type dr  = 1e-5;

        test::check_force(sys, forcefield, tol, dr);
        test::check_virial(sys, forcefield, tol);
    }
}

BOOST_AUTO_TEST_CASE(MultipleBasin_2Basin_numerical_difference_bias)
{
    namespace test = mjolnir::test;

    mjolnir::LoggerManager::set_default_logger("test_multiple_basin_forcefield.log");

    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
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
            /*k_chi = */1.0, /*chi_0 = */ 0.5);

    std::vector<std::unique_ptr<mjolnir::MultipleBasinUnitBase<traits_type>>> units;
    units.push_back(std::move(unit1));

    std::mt19937 mt(123456789);

    mjolnir::MultipleBasinForceField<traits_type> forcefield(
            std::make_tuple(std::move(locc), std::move(gloc), std::move(extc)),
            std::move(units));

    for(const real_type len : {1.5, 1.75, 2.0})
    {
        for(std::size_t i = 0; i < 1000; ++i)
        {
            system_type sys(2, boundary_type{});

            test::clear_everything(sys);

            sys.position(0) = coord_type( 0.0, 0.0, 0.0);
            sys.position(1) = coord_type( len, 0.0, 0.0);

            test::apply_random_rotation(sys, mt);
            test::apply_random_perturbation(sys, mt, 0.01);

            forcefield.initialize(sys); // don't forget this

            constexpr real_type tol = 1e-4;
            constexpr real_type dr  = 1e-5;

            // Here we need to take sum of `chi_bias` and `V_MB`.
            // So we don't use test::check_force.

            const auto init = sys; // save current system state as the initial state
            for(std::size_t idx=0; idx<sys.size(); ++idx)
            {
                {
                    // ----------------------------------------------------------------
                    // reset positions
                    sys = init;

                    // calc U(x-dx) and U(x+dx)
                    auto E0 = forcefield.calc_force_and_energy(sys);
                    E0 += sys.attribute("chi_bias");

                    // calc U(x+dx)
                    mjolnir::math::X(sys.position(idx)) += 2 * dr;
                    auto E1 = forcefield.calc_force_and_energy(sys);
                    E1 += sys.attribute("chi_bias");

                    // take the center
                    sys = init;
                    mjolnir::math::X(sys.position(idx)) += dr;

                    // calc F(x)
                    sys.preprocess_forces();
                    forcefield.calc_force(sys);
                    sys.postprocess_forces();

                    // central difference
                    const auto dE = (E1 - E0) * 0.5;

                    BOOST_TEST(-dE == dr * mjolnir::math::X(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                }
                {
                    // ----------------------------------------------------------------
                    // reset positions
                    sys = init;

                    // calc U(x-dx) and U(x+dx)
                    auto E0 = forcefield.calc_force_and_energy(sys);
                    E0 += sys.attribute("chi_bias");

                    // calc U(x+dx)
                    mjolnir::math::Y(sys.position(idx)) += 2 * dr;
                    auto E1 = forcefield.calc_force_and_energy(sys);
                    E1 += sys.attribute("chi_bias");

                    // take the center
                    sys = init;
                    mjolnir::math::Y(sys.position(idx)) += dr;

                    // calc F(x)
                    sys.preprocess_forces();
                    forcefield.calc_force(sys);
                    sys.postprocess_forces();

                    // central difference
                    const auto dE = (E1 - E0) * 0.5;

                    BOOST_TEST(-dE == dr * mjolnir::math::Y(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                }
                {
                    // ----------------------------------------------------------------
                    // reset positions
                    sys = init;

                    // calc U(x-dx) and U(x+dx)
                    auto E0 = forcefield.calc_force_and_energy(sys);
                    E0 += sys.attribute("chi_bias");

                    // calc U(x+dx)
                    mjolnir::math::Z(sys.position(idx)) += 2 * dr;
                    auto E1 = forcefield.calc_force_and_energy(sys);
                    E1 += sys.attribute("chi_bias");

                    // take the center
                    sys = init;
                    mjolnir::math::Z(sys.position(idx)) += dr;

                    // calc F(x)
                    sys.preprocess_forces();
                    forcefield.calc_force(sys);
                    sys.postprocess_forces();

                    // central difference
                    const auto dE = (E1 - E0) * 0.5;

                    BOOST_TEST(-dE == dr * mjolnir::math::Z(sys.force(idx)),
                               boost::test_tools::tolerance(tol));
                }
            }
        }
    }
}


BOOST_AUTO_TEST_CASE(MultipleBasin_3Basin_numerical_difference)
{
    namespace test = mjolnir::test;

    mjolnir::LoggerManager::set_default_logger("test_multiple_basin_forcefield.log");

    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
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

    for(std::size_t i = 0; i < 1000; ++i)
    {
        system_type sys(2, boundary_type{});

        test::clear_everything(sys);

        sys.position(0) = coord_type( 0.0, 0.0, 0.0);
        sys.position(1) = coord_type( 1.0, 1.0, 1.0);

        test::apply_random_rotation(sys, mt);
        test::apply_random_perturbation(sys, mt, 0.01);

        forcefield.initialize(sys); // don't forget this

        constexpr real_type tol = 1e-4;
        constexpr real_type dr  = 1e-5;
        test::check_force(sys, forcefield, tol, dr);
        test::check_virial(sys, forcefield, tol);
    }
}
