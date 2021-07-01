#define BOOST_TEST_MODULE "test_excluded_volume_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>
#include <mjolnir/forcefield/3SPN2/ThreeSPN2ExcludedVolumePotential.hpp>

BOOST_AUTO_TEST_CASE(EXV_double)
{
    mjolnir::LoggerManager::set_default_logger("test_3spn2_excluded_volume_potential.log");
    using real_type = double;
    using potential_type = mjolnir::ThreeSPN2ExcludedVolumePotential<real_type>;
    using parameter_type = potential_type::parameter_type;
    using traits_type    = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using boundary_type  = typename traits_type::boundary_type;
    using system_type    = mjolnir::System<traits_type>;

    {
        // to make the value realistic, we need to modify the unit system.
        using phys_type = mjolnir::physics::constants<real_type>;
        using unit_type = mjolnir::unit::constants<real_type>;

        phys_type::set_energy_unit("kcal/mol");
        phys_type::set_length_unit("angstrom");

        phys_type::set_kB(unit_type::boltzmann_constant() *
                          1e-3 * unit_type::J_to_cal() * unit_type::avogadro_constant());
        phys_type::set_NA(unit_type::avogadro_constant());
        phys_type::set_eps0((unit_type::vacuum_permittivity() /
            unit_type::elementary_charge()) / unit_type::elementary_charge() *
            (1e+3 / unit_type::J_to_cal() / unit_type::avogadro_constant()) *
            (1.0 / unit_type::m_to_angstrom()));

        phys_type::set_m_to_length(1e-10);
        phys_type::set_length_to_m(1e+10);
        phys_type::set_L_to_volume(1e+27);
        phys_type::set_volume_to_L(1e-27);
    }

    constexpr std::size_t N = 10000;
    constexpr real_type   h = 1e-6;
    constexpr real_type tol = 1e-6;

    const real_type sigma   = 3.0;

    potential_type potential;

    system_type sys(1, boundary_type{});
    potential.initialize(sys);

    const real_type x_min = 0.8 * sigma;
    const real_type x_max = potential.cutoff_ratio() * sigma;

    mjolnir::test::check_potential(potential, parameter_type{sigma},
                                   x_min, x_max, tol, h, N);
}

BOOST_AUTO_TEST_CASE(EXV_float)
{
    mjolnir::LoggerManager::set_default_logger("test_3spn2_excluded_volume_potential.log");
    using real_type = float;
    using potential_type = mjolnir::ThreeSPN2ExcludedVolumePotential<real_type>;
    using parameter_type = potential_type::parameter_type;
    using traits_type    = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using boundary_type  = typename traits_type::boundary_type;
    using system_type    = mjolnir::System<traits_type>;

    {
        // to make the value realistic, we need to modify the unit system.
        using phys_type = mjolnir::physics::constants<real_type>;
        using unit_type = mjolnir::unit::constants<real_type>;

        phys_type::set_energy_unit("kcal/mol");
        phys_type::set_length_unit("angstrom");

        phys_type::set_kB(unit_type::boltzmann_constant() *
                          1e-3 * unit_type::J_to_cal() * unit_type::avogadro_constant());
        phys_type::set_NA(unit_type::avogadro_constant());
        phys_type::set_eps0((unit_type::vacuum_permittivity() /
            unit_type::elementary_charge()) / unit_type::elementary_charge() *
            (1e+3 / unit_type::J_to_cal() / unit_type::avogadro_constant()) *
            (1.0 / unit_type::m_to_angstrom()));

        phys_type::set_m_to_length(1e-10);
        phys_type::set_length_to_m(1e+10);
        phys_type::set_L_to_volume(1e+27);
        phys_type::set_volume_to_L(1e-27);
    }

    constexpr static std::size_t N = 200;
    constexpr static real_type   h = 0.002f;
    constexpr static real_type tol = 0.005f;

    const real_type sigma   = 3.0f;

    potential_type potential;

    system_type sys(1, boundary_type{});
    potential.initialize(sys);

    const real_type x_min = 0.8f * sigma;
    const real_type x_max = potential.cutoff_ratio() * sigma;

    mjolnir::test::check_potential(potential, parameter_type{sigma},
                                   x_min, x_max, tol, h, N);
}

BOOST_AUTO_TEST_CASE(EXV_cutoff_length)
{
    mjolnir::LoggerManager::set_default_logger("test_3spn2_excluded_volume_potential.log");
    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using boundary_type  = typename traits_type::boundary_type;
    using system_type    = mjolnir::System<traits_type>;

    using potential_type = mjolnir::ThreeSPN2ExcludedVolumePotential<real_type>;
    using parameter_list_type  = mjolnir::ThreeSPN2ExcludedVolumeParameterList<traits_type>;
    using parameter_type       = parameter_list_type::parameter_type;
    using ignore_molecule_type = parameter_list_type::ignore_molecule_type;
    using ignore_group_type    = parameter_list_type::ignore_group_type;

    {
        // to make the value realistic, we need to modify the unit system.
        using phys_type = mjolnir::physics::constants<real_type>;
        using unit_type = mjolnir::unit::constants<real_type>;

        phys_type::set_energy_unit("kcal/mol");
        phys_type::set_length_unit("angstrom");

        phys_type::set_kB(unit_type::boltzmann_constant() *
                          1e-3 * unit_type::J_to_cal() * unit_type::avogadro_constant());
        phys_type::set_NA(unit_type::avogadro_constant());
        phys_type::set_eps0((unit_type::vacuum_permittivity() /
            unit_type::elementary_charge()) / unit_type::elementary_charge() *
            (1e+3 / unit_type::J_to_cal() / unit_type::avogadro_constant()) *
            (1.0 / unit_type::m_to_angstrom()));

        phys_type::set_m_to_length(1e-10);
        phys_type::set_length_to_m(1e+10);
        phys_type::set_L_to_volume(1e+27);
        phys_type::set_volume_to_L(1e-27);
    }

    potential_type potential;
    system_type sys(1, boundary_type{});
    potential.initialize(sys);

    // generate a parameter
    std::vector<std::pair<std::size_t, parameter_type>> parameters(10);
    for(std::size_t i=0; i<parameters.size(); ++i)
    {
        parameters.at(i).first = i;
    }
    parameters.at(0).second = mjolnir::parameter_3SPN2::bead_kind::BaseA;
    parameters.at(1).second = mjolnir::parameter_3SPN2::bead_kind::BaseT;
    parameters.at(2).second = mjolnir::parameter_3SPN2::bead_kind::BaseC;
    parameters.at(3).second = mjolnir::parameter_3SPN2::bead_kind::BaseG;
    parameters.at(4).second = mjolnir::parameter_3SPN2::bead_kind::BaseA;
    parameters.at(5).second = mjolnir::parameter_3SPN2::bead_kind::BaseT;
    parameters.at(6).second = mjolnir::parameter_3SPN2::bead_kind::BaseC;
    parameters.at(7).second = mjolnir::parameter_3SPN2::bead_kind::BaseG;
    parameters.at(8).second = mjolnir::parameter_3SPN2::bead_kind::Phosphate;
    parameters.at(9).second = mjolnir::parameter_3SPN2::bead_kind::Sugar;

    parameter_list_type parameter_list(
        mjolnir::ThreeSPN2ExcludedVolumePotentialParameter<real_type>{},
        parameters, {}, ignore_molecule_type{"Nothing"}, ignore_group_type({}));

    // secure maximum cutoff calculated from per-pair parameter
    // This requires O(N^2), so it will not calculated in usual cases
    real_type max_abs_cutoff = 0.0;
    for(std::size_t i=0; i<parameters.size(); ++i)
    {
        for(std::size_t j=0; j<parameters.size(); ++j)
        {
            max_abs_cutoff = std::max(max_abs_cutoff, potential.absolute_cutoff(parameter_list.prepare_params(i, j)));
        }
    }

    const real_type safety = 1.0 + std::numeric_limits<real_type>::epsilon();
    for(std::size_t i=0; i<parameters.size(); ++i)
    {
        for(std::size_t j=0; j<parameters.size(); ++j)
        {
            BOOST_TEST(potential.potential(max_abs_cutoff * safety, parameter_list.prepare_params(i, j)) == 0.0);
        }
    }
}
