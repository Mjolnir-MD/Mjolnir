#define BOOST_TEST_MODULE "test_debye_huckel_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/utility.hpp>

#include <mjolnir/forcefield/global/DebyeHuckelPotential.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>

BOOST_AUTO_TEST_CASE(DH_double)
{
    using real_type = double;
    using potential_type = mjolnir::DebyeHuckelPotential<real_type>;
    using parameter_type = potential_type::parameter_type;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using boundary_type = typename traits_type::boundary_type;

    mjolnir::LoggerManager::set_default_logger("test_debye_huckel_potential.log");
    {
        // to make the value realistic, we need to modify the unit system.
        using phys_type = mjolnir::physics::constants<real_type>;
        using unit_type = mjolnir::unit::constants<real_type>;

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

    mjolnir::System<traits_type> sys(0, boundary_type{});
    sys.attribute("temperature")    = 300.0; // [K]
    sys.attribute("ionic_strength") =   0.2; // [M]

    potential_type potential(5.5);
    potential.initialize(sys);

    parameter_type params{1.0};

    const real_type x_min = 0.5 * potential.debye_length();
    const real_type x_max = 5.5 * potential.debye_length();

    mjolnir::test::check_potential(potential, params, x_min, x_max, tol, h, N);
}
BOOST_AUTO_TEST_CASE(DH_float)
{
    using real_type = float;
    using potential_type = mjolnir::DebyeHuckelPotential<real_type>;
    using parameter_type = potential_type::parameter_type;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using boundary_type = typename traits_type::boundary_type;

    mjolnir::LoggerManager::set_default_logger("test_debye_huckel_potential.log");
    {
        // to make the value realistic, we need to modify the unit system.
        using phys_type = mjolnir::physics::constants<real_type>;
        using unit_type = mjolnir::unit::constants<real_type>;

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

    constexpr static std::size_t N = 1000;
    constexpr static real_type   h = 0.002f;
    constexpr static real_type tol = 0.004f;

    mjolnir::System<traits_type> sys(0, boundary_type{});
    sys.attribute("temperature")    = 300.0; // [K]
    sys.attribute("ionic_strength") =   0.2; // [M]

    potential_type potential(5.5f);
    potential.initialize(sys);

    parameter_type params{1.0f};

    const real_type x_min = 0.5f * potential.debye_length();
    const real_type x_max = 5.5f * potential.debye_length();

    mjolnir::test::check_potential(potential, params, x_min, x_max, tol, h, N);
}

BOOST_AUTO_TEST_CASE(DH_cutoff_length)
{
    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using boundary_type = typename traits_type::boundary_type;
    using potential_type = mjolnir::DebyeHuckelPotential<real_type>;

    using parameter_list_type  = mjolnir::DebyeHuckelParameterList<traits_type>;
    using parameter_type       = parameter_list_type::parameter_type;
    using ignore_molecule_type = parameter_list_type::ignore_molecule_type;
    using ignore_group_type    = parameter_list_type::ignore_group_type;

    mjolnir::LoggerManager::set_default_logger("test_debye_huckel_potential.log");
    {
        // to make the value realistic, we need to modify the unit system.
        using phys_type = mjolnir::physics::constants<real_type>;
        using unit_type = mjolnir::unit::constants<real_type>;

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

    mjolnir::System<traits_type> sys(10, boundary_type{});
    sys.attribute("temperature")    = 300.0; // [K]
    sys.attribute("ionic_strength") =   0.2; // [M]

    potential_type potential(5.5f);
    potential.initialize(sys);

    // generate a parameter
    std::vector<std::pair<std::size_t, parameter_type>> parameters(10);
    for(std::size_t i=0; i<parameters.size(); ++i)
    {
        parameters.at(i).first          = i;
        parameters.at(i).second.charge  = (i % 2 == 0) ? 1.0 : -1.0;
    }
    parameter_list_type parameter_list(parameters, {}, ignore_molecule_type{"Nothing"}, ignore_group_type({}));

    // estimated maximum cutoff (from per-particle parameter)
    const auto max_cutoff = potential.max_cutoff(parameter_list.parameters().begin(), parameter_list.parameters().end());

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

    // max cutoff estimated from per-particle parameter should be the same or larger than the secure max cutoff
    BOOST_TEST(max_abs_cutoff <= max_cutoff);

    const real_type safety = 1.0 + std::numeric_limits<real_type>::epsilon();
    for(std::size_t i=0; i<parameters.size(); ++i)
    {
        for(std::size_t j=0; j<parameters.size(); ++j)
        {
            BOOST_TEST(potential.potential(max_abs_cutoff * safety, parameter_list.prepare_params(i, j)) == 0.0);
        }
    }
}
