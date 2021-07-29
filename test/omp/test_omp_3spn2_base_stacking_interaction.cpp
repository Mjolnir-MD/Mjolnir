#define BOOST_TEST_MODULE "test_omp_3spn2_base_stacking_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/math/constants.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/RandomNumberGenerator.hpp>
#include <mjolnir/omp/ThreeSPN2BaseStackingInteraction.hpp>
#include <mjolnir/omp/ThreeSPN2BaseStackingInteraction.hpp>
#include <mjolnir/input/read_units.hpp>

BOOST_AUTO_TEST_CASE(OpenMP_ThreeSPN2BaseStackingInteraction)
{
    namespace test = mjolnir::test;

    mjolnir::LoggerManager::set_default_logger(
            "test_omp_3spn2_base_stacking_interaction.log");

    using seq_traits_type = mjolnir::SimulatorTraits      <double, mjolnir::UnlimitedBoundary>;
    using omp_traits_type = mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>;

    using real_type       = typename omp_traits_type::real_type;
    using coordinate_type = typename omp_traits_type::coordinate_type;
    using boundary_type   = typename omp_traits_type::boundary_type;

    using seq_system_type = mjolnir::System<seq_traits_type>;
    using omp_system_type = mjolnir::System<omp_traits_type>;

    using potential_type  = mjolnir::ThreeSPN2BaseStackingPotential<real_type>;
    using base_stack_kind = typename potential_type::base_stack_kind;

    using omp_interaction_type = mjolnir::ThreeSPN2BaseStackingInteraction<omp_traits_type>;
    using seq_interaction_type = mjolnir::ThreeSPN2BaseStackingInteraction<seq_traits_type>;

    using ParameterSet = mjolnir::ThreeSPN2BaseStackingPotentialParameter<real_type>;
    {
        using unit_type = mjolnir::unit::constants<real_type>;
        using phys_type = mjolnir::physics::constants<real_type>;
        const std::string energy = "kcal/mol";
        const std::string length = "angstrom";
        phys_type::set_kB(phys_type::kB() * (unit_type::J_to_cal() / 1000.0) *
                          unit_type::avogadro_constant());
        phys_type::set_eps0(phys_type::eps0() * (1000.0 / unit_type::J_to_cal()) /
                            unit_type::avogadro_constant());
        phys_type::set_energy_unit(energy);

        phys_type::set_eps0(phys_type::eps0() / unit_type::m_to_angstrom());

        phys_type::set_m_to_length(unit_type::m_to_angstrom());
        phys_type::set_length_to_m(unit_type::angstrom_to_m());

        phys_type::set_L_to_volume(1e-3 * std::pow(unit_type::m_to_angstrom(), 3));
        phys_type::set_volume_to_L(1e+3 * std::pow(unit_type::angstrom_to_m(), 3));

        phys_type::set_length_unit(length);
    }

    //        SBi
    //     Si --> Bi
    //    /     `-^
    //   Pj theta | rij
    //    \       |
    //     Sj --- Bj
    //
    //  rij:
    //  1. (r < r0)
    //  2. (r0 < r)
    //  theta:
    //  1. theta < pi/2K
    //  2. pi/2K < theta < pi/K
    //  3. pi/K  < theta

    std::mt19937 rng(123456789);

    const auto bs_kind = base_stack_kind::AA;

    potential_type   potential(ParameterSet{});
    omp_interaction_type omp_interaction("none",
            std::vector<std::pair<std::array<std::size_t, 3>, base_stack_kind>>{
                { {{0, 1, 2}}, bs_kind }
            }, potential_type(potential), {});
    seq_interaction_type seq_interaction("none",
            std::vector<std::pair<std::array<std::size_t, 3>, base_stack_kind>>{
                { {{0, 1, 2}}, bs_kind }
            }, potential_type(potential), {});

    omp_system_type omp_sys(3, boundary_type{});
    seq_system_type seq_sys(3, boundary_type{});

    test::clear_everything(omp_sys);
    test::clear_everything(seq_sys);

    omp_sys.name(0)  = "Si";
    omp_sys.name(1)  = "Bi";
    omp_sys.name(2)  = "Bj";
    omp_sys.group(0) = "DNA";
    omp_sys.group(1) = "DNA";
    omp_sys.group(2) = "DNA";

    potential.initialize(omp_sys);

    omp_interaction.initialize(omp_sys);
    seq_interaction.initialize(seq_sys);

    const auto theta0    = potential.theta_0(bs_kind);
    const auto pi_over_K = potential.pi_over_K_BS();
    const auto theta0_1  = theta0 + 0.2 * pi_over_K; //         dtheta < pi/2K
    const auto theta0_2  = theta0 + 0.7 * pi_over_K; // pi/2K < dtheta < pi/K
    const auto theta0_3  = theta0 + 1.2 * pi_over_K; // pi/K  < dtheta
    const auto r0_1      = potential.r0(bs_kind) - 0.2;
    const auto r0_2      = potential.r0(bs_kind) + 0.5;

    for(const auto r : {r0_1, r0_2})
    {
    for(const auto theta : {theta0_1, theta0_2, theta0_3})
    {
        BOOST_TEST_MESSAGE("======================================");
        BOOST_TEST_MESSAGE("r = " << r << ", theta = " << theta);

    for(std::size_t i=0; i<100; ++i)
    {
        // generate particle configuration in the following way
        //    y
        // Bj ^
        //  \ | theta0
        // r0\|-.
        // ---o-----o--> x
        //  Bi      Si
        //

        omp_sys.position(0) = coordinate_type(4.0, 0.0, 0.0); // Si
        omp_sys.position(1) = coordinate_type(0.0, 0.0, 0.0); // Bi
        omp_sys.position(2) = coordinate_type(r * std::cos(theta), r * std::sin(theta), 0.0); // Bj

        test::apply_random_rotation(omp_sys, rng);
        test::apply_random_perturbation(omp_sys, rng, 0.01);
        test::clear_force(omp_sys);
        test::clear_force(seq_sys);

        seq_sys.position(0) = omp_sys.position(0);
        seq_sys.position(1) = omp_sys.position(1);
        seq_sys.position(2) = omp_sys.position(2);

        constexpr real_type tol = 1e-4;

        test::check_force_consistency              (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
        test::check_force_and_energy_consistency   (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
        test::check_force_and_virial_consistency   (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
        test::check_force_energy_virial_consistency(omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
        test::check_energy_consistency             (omp_sys, omp_interaction, seq_sys, seq_interaction, tol);
    } // perturbation
    } // theta
    } // r
}
