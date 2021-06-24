#define BOOST_TEST_MODULE "test_lennard_jones_attractive_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/util/check_potential.hpp>
#include <mjolnir/forcefield/global/ParameterList.hpp>
#include <mjolnir/forcefield/global/LennardJonesPotential.hpp>
#include <mjolnir/forcefield/global/LennardJonesAttractivePotential.hpp>

BOOST_AUTO_TEST_CASE(LennardJones_double)
{
    using real_type = double;
    using potential_type = mjolnir::LennardJonesAttractivePotential<real_type>;
    using parameter_type = potential_type::parameter_type;

    constexpr std::size_t N = 10000;
    constexpr real_type   h = 1e-6;
    constexpr real_type tol = 1e-6;

    const real_type sigma   = 3.0;
    const real_type epsilon = 1.0;
    const parameter_type para{sigma, epsilon};

    potential_type potential(/*cutoff = */ 2.5);

    const real_type x_min = 0.8 * sigma;
    const real_type x_max = potential.cutoff_ratio() * sigma;

    mjolnir::test::check_potential(potential, para, x_min, x_max, tol, h, N);
    for(std::size_t i=1; i<N; ++i)
    {
        const real_type x = (0.8 + (std::pow(2.0, 1.0/6.0) - 0.8) * i / N) * sigma;
        BOOST_TEST(potential.potential(x, para)  == -epsilon, boost::test_tools::tolerance(h));
        BOOST_TEST(potential.derivative(x, para) == 0.0, boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(LennardJones_float)
{
    using real_type = float;
    using potential_type = mjolnir::LennardJonesAttractivePotential<real_type>;
    using parameter_type = potential_type::parameter_type;

    constexpr std::size_t N = 1000;
    constexpr real_type   h = 0.002f;
    constexpr real_type tol = 0.005f;

    const real_type sigma   = 3.0f;
    const real_type epsilon = 1.0f;
    const parameter_type para{sigma, epsilon};

    potential_type potential(/*cutoff = */ 2.5f);

    const real_type x_min = std::pow(2.0f, 1.0f/6.0f) * sigma;
    const real_type x_max = potential.cutoff_ratio() * sigma;

    mjolnir::test::check_potential(potential, parameter_type{sigma, epsilon},
                                   x_min, x_max, tol, h, N);

    for(std::size_t i=1; i<N; ++i)
    {
        const real_type x = (0.8f + (std::pow(2.0f, 1.0f/6.0f) - 0.8f) * i / N) * sigma;
        BOOST_TEST(potential.potential(x, para)  == -epsilon, boost::test_tools::tolerance(h));
        BOOST_TEST(potential.derivative(x, para) == 0.0f, boost::test_tools::tolerance(h));
    }
}

BOOST_AUTO_TEST_CASE(LJ_cutoff_length)
{
    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::LennardJonesAttractivePotential<real_type>;

    using parameter_list_type  = mjolnir::LorentzBerthelotRule<traits_type, potential_type>;
    using parameter_type       = parameter_list_type::parameter_type;
    using ignore_molecule_type = parameter_list_type::ignore_molecule_type;
    using ignore_group_type    = parameter_list_type::ignore_group_type;

    potential_type potential(2.5);

    // generate a parameter
    std::vector<std::pair<std::size_t, parameter_type>> parameters(10);
    for(std::size_t i=0; i<parameters.size(); ++i)
    {
        parameters.at(i).first          = i;
        parameters.at(i).second.epsilon = 1.0;
        parameters.at(i).second.sigma   = static_cast<real_type>(i) + 1.0;
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
