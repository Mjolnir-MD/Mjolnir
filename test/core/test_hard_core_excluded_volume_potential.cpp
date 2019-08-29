#define BOOST_TEST_MODULE "test_excluded_volume_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/potential/global/HardCoreExcludedVolumePotential.hpp>

BOOST_AUTO_TEST_CASE(HCEXV_double)
{
    using real_type = double;
    using molecule_id_type = mjolnir::Topology::molecule_id_type;
    using group_id_type    = mjolnir::Topology::group_id_type;
    using potential_type   = mjolnir::HardCoreExcludedVolumePotential<real_type>;
    using parameter_type   = potential_type::parameter_type;

    constexpr static std::size_t N = 10000;
    constexpr static real_type   h = 1e-6;

    const real_type epsilon = 1.0;
    const real_type hard_core_radius = 10.0;
    const real_type soft_shell_thickness = 3.0;
    const parameter_type param{hard_core_radius, soft_shell_thickness};
    mjolnir::HardCoreExcludedVolumePotential<real_type> hdexv{
      epsilon, {{0, param}, {1, param}}, {},
      mjolnir::IgnoreMolecule<molecule_id_type>("Nothing"),
      mjolnir::IgnoreGroup<group_id_type>({})
    };

    const real_type x_min = 0.8 * soft_shell_thickness + hard_core_radius * 2;
    const real_type x_max = hdexv.cutoff_ratio() * soft_shell_thickness + hard_core_radius * 2;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
      const real_type x = x_min + i * dx;
      const real_type pot1 = hdexv.potential(0, 1, x + h);
      const real_type pot2 = hdexv.potential(0, 1, x - h);
      const real_type dpot = (pot1 - pot2) / (2 * h);
      const real_type deri = hdexv.derivative(0, 1, x);

      BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(HCEXV_float)
{
    using real_type = float;
    using molecule_id_type = mjolnir::Topology::molecule_id_type;
    using group_id_type    = mjolnir::Topology::group_id_type;
    using potential_type   = mjolnir::HardCoreExcludedVolumePotential<real_type>;
    using parameter_type   = potential_type::parameter_type;

    constexpr static std::size_t N = 1000;
    constexpr static real_type   h = 0.002f;
    constexpr static real_type tol = 0.005f;

    const real_type epsilon = 1.0;
    const real_type hard_core_radius = 10.0;
    const real_type soft_shell_thickness = 3.0;
    const parameter_type param{hard_core_radius, soft_shell_thickness};
    mjolnir::HardCoreExcludedVolumePotential<real_type> hdexv{
      epsilon, {{0, param}, {1, param}}, {},
      mjolnir::IgnoreMolecule<molecule_id_type>("Nothing"),
      mjolnir::IgnoreGroup<group_id_type>({})
    };

    const real_type x_min = 0.8 * soft_shell_thickness + hard_core_radius * 2;
    const real_type x_max = hdexv.cutoff_ratio() * soft_shell_thickness + hard_core_radius * 2;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
      const real_type x = x_min + i * dx;
      const real_type pot1 = hdexv.potential(0, 1, x + h);
      const real_type pot2 = hdexv.potential(0, 1, x - h);
      const real_type dpot = (pot1 - pot2) / (2 * h);
      const real_type deri = hdexv.derivative(0, 1, x);

      BOOST_TEST(dpot == deri, boost::test_tools::tolerance(tol));
    }
}
