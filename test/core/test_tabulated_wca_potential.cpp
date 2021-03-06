#define BOOST_TEST_MODULE "test_tabulated_wca_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/forcefield/global/TabulatedWCAPotential.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

BOOST_AUTO_TEST_CASE(WCA_double)
{
    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = typename traits_type::real_type;
    using molecule_id_type = mjolnir::Topology::molecule_id_type;
    using group_id_type    = mjolnir::Topology::group_id_type;
    using potential_type   = mjolnir::TabulatedWCAPotential<traits_type>;
    using parameter_type   = potential_type::pair_parameter_type;

    constexpr static std::size_t N = 10000;
    constexpr static real_type   h = 1e-6;

    const real_type sigma   = 3.0;
    const real_type epsilon = 1.0;
    const parameter_type param{sigma, epsilon};
    potential_type wca{
        potential_type::default_cutoff(),
        {{"A:B", {sigma, epsilon}}},
        {{0, "A"}, {1, "B"}}, {},
        mjolnir::IgnoreMolecule<molecule_id_type>("Nothing"),
        mjolnir::IgnoreGroup   <group_id_type   >({})
    };

    const real_type x_min = 0.8 * sigma;
    const real_type x_max = wca.cutoff_ratio() * sigma - 2 * h;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = wca.potential(0, 1, x + h);
        const real_type pot2 = wca.potential(0, 1, x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = wca.derivative(0, 1, x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
    }

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x   = x_max + i * dx;
        const real_type pot = wca.potential(0, 1, x);
        BOOST_TEST(pot == 0.0, boost::test_tools::tolerance(h));
    }
}
