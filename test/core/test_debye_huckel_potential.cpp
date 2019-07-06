#define BOOST_TEST_MODULE "test_debye_huckel_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/potential/global/DebyeHuckelPotential.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/logger.hpp>

BOOST_AUTO_TEST_CASE(DH_double)
{
    mjolnir::LoggerManager::set_default_logger("test_debye_huckel_potential.log");

    using real_type = double;
    using molecule_id_type = mjolnir::Topology::molecule_id_type;
    using group_id_type    = mjolnir::Topology::group_id_type;

    // to make the value realistic, we need to modify the unit system.
    using phys_type = mjolnir::physics::constants<real_type>;
    using unit_type = mjolnir::unit::constants<real_type>;

    phys_type::set_kB(unit_type::boltzmann_constant *
                      1e-3 * unit_type::J_to_cal * unit_type::avogadro_constant);
    phys_type::set_NA(unit_type::avogadro_constant);
    phys_type::set_eps0((unit_type::vacuum_permittivity /
        unit_type::elementary_charge) / unit_type::elementary_charge *
        (1e+3 / unit_type::J_to_cal / unit_type::avogadro_constant) *
        (1.0 / unit_type::m_to_angstrom));

    phys_type::set_m_to_length(1e-10);
    phys_type::set_length_to_m(1e+10);
    phys_type::set_L_to_volume(1e+27);
    phys_type::set_volume_to_L(1e-27);

    constexpr std::size_t N = 10000;
    constexpr real_type   h = 1e-6;

    const real_type charge = 1.0;
    mjolnir::DebyeHuckelPotential<real_type> dh(
        {{0u, charge}, {1u, charge}}, {},
        mjolnir::IgnoreMolecule<molecule_id_type>("Nothing"),
        mjolnir::IgnoreGroup   <group_id_type   >("Nothing")
        );

    const real_type x_min = 0.5 * dh.debye_length();
    const real_type x_max = dh.max_cutoff_length();
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = dh.potential(0, 1, x + h);
        const real_type pot2 = dh.potential(0, 1, x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = dh.derivative(0, 1, x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(DH_float)
{
    using real_type = float;
    using molecule_id_type = mjolnir::Topology::molecule_id_type;
    using group_id_type    = mjolnir::Topology::group_id_type;

    using phys_type = mjolnir::physics::constants<real_type>;
    using unit_type = mjolnir::unit::constants<real_type>;

    phys_type::set_kB(unit_type::boltzmann_constant *
        (1e-3 * unit_type::J_to_cal * unit_type::avogadro_constant));
    phys_type::set_NA(unit_type::avogadro_constant);
    phys_type::set_eps0((unit_type::vacuum_permittivity /
        unit_type::elementary_charge) / unit_type::elementary_charge *
        (1e+3 / unit_type::J_to_cal / unit_type::avogadro_constant) *
        (1.0 / unit_type::m_to_angstrom));

    phys_type::set_m_to_length(1e-10);
    phys_type::set_length_to_m(1e+10);
    phys_type::set_L_to_volume(1e+27);
    phys_type::set_volume_to_L(1e-27);

    constexpr static std::size_t N = 1000;
    constexpr static real_type   h = 1e-2;

    const real_type charge = 1.0;
    mjolnir::DebyeHuckelPotential<real_type> dh(
        {{0u, charge}, {1u, charge}}, {},
        mjolnir::IgnoreMolecule<molecule_id_type>("Nothing"),
        mjolnir::IgnoreGroup   <group_id_type   >("Nothing")
        );

    const real_type x_min = 0.5 * dh.debye_length();
    const real_type x_max = dh.max_cutoff_length();
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = dh.potential(0, 1, x + h);
        const real_type pot2 = dh.potential(0, 1, x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = dh.derivative(0, 1, x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
    }
}
