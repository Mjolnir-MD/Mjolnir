#define BOOST_TEST_MODULE "test_isolf_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/forcefield/iSoLF/iSoLFAttractivePotential.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

BOOST_AUTO_TEST_CASE(iSoLF_double)
{
    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type   = typename traits_type::real_type;
    using molecule_id_type = mjolnir::Topology::molecule_id_type;
    using group_id_type    = mjolnir::Topology::group_id_type;
    using potential_type   = mjolnir::iSoLFAttractivePotential<traits_type>;
    using parameter_type   = potential_type::parameter_type;

    constexpr static std::size_t N = 10000;
    constexpr static real_type   h = 1e-6;

    const real_type sigma   = 3.0;
    const real_type epsilon = 1.0;
    const real_type omega   = 1.0;
    const parameter_type param{sigma, epsilon, omega};
    potential_type isolf{
        {{0, param}, {1, param}}, {},
        mjolnir::IgnoreMolecule<molecule_id_type>("Nothing"),
        mjolnir::IgnoreGroup   <group_id_type   >({})
    };

    const real_type x_min = 1.3 * sigma;
    const real_type x_max = isolf.max_cutoff_length() * 0.9;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = isolf.potential(0, 1, x + h);
        const real_type pot2 = isolf.potential(0, 1, x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = isolf.derivative(0, 1, x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(iSoLF_float)
{
    using traits_type = mjolnir::SimulatorTraits<float, mjolnir::UnlimitedBoundary>;
    using real_type   = typename traits_type::real_type;
    using molecule_id_type = mjolnir::Topology::molecule_id_type;
    using group_id_type    = mjolnir::Topology::group_id_type;
    using potential_type   = mjolnir::iSoLFAttractivePotential<traits_type>;
    using parameter_type   = potential_type::parameter_type;
    constexpr std::size_t N = 1000;
    constexpr real_type   h = 0.002f;
    constexpr real_type tol = 0.005f;

    const real_type sigma   = 3.0f;
    const real_type epsilon = 1.0f;
    const real_type omega   = 1.0f;
    const parameter_type param{sigma, epsilon, omega};
    potential_type isolf{
        {{0, param}, {1, param}}, {},
        mjolnir::IgnoreMolecule<molecule_id_type>("Nothing"),
        mjolnir::IgnoreGroup   <group_id_type   >({})
    };

    const real_type x_min = 1.3f * sigma;
    const real_type x_max = isolf.max_cutoff_length() * 0.9;
    const real_type dx = (x_max - x_min) / N;

    for(std::size_t i=0; i<N; ++i)
    {
        const real_type x    = x_min + i * dx;
        const real_type pot1 = isolf.potential(0, 1, x + h);
        const real_type pot2 = isolf.potential(0, 1, x - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = isolf.derivative(0, 1, x);

        BOOST_TEST(dpot == deri, boost::test_tools::tolerance(tol));
    }
}


