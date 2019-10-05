#define BOOST_TEST_MODULE "test_read_directional_contact_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/input/read_local_interaction.hpp>

BOOST_AUTO_TEST_CASE(read_directional_contact_cosine_go_contact)
{
    mjolnir::LoggerManager::set_default_logger("test_read_directional_contact_interaction.log");
    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type = traits_type::real_type;
    const real_type tol = 1e-7;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
          interaction        = "DirectionalContact"
          potentials.angle1   = "Cosine"
          potentials.angle2   = "Cosine"
          potentials.contact  = "GoContact"
          topology           = "none"
          parameters = [
          {indices = [0, 1, 2, 3], angle1 = {v0 = 0.0, k = -10.0, n = 1}, angle2 = {v0 = 1.0, k = -20.0, n = 2}, contact = {v0 = 2.0, k = 1.0}},
          ]
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv =
                dynamic_cast<
                    mjolnir::DirectionalContactInteraction<
                        traits_type, mjolnir::CosinePotential<real_type>,
                        mjolnir::CosinePotential<real_type>,
                        mjolnir::GoContactPotential<real_type>
                        >*
                >(base.get()); //check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
        const auto indices_potentials_tuple = derv->potentials().at(0);
        const auto angle1_pot  = std::get<1>(indices_potentials_tuple);
        const auto angle2_pot  = std::get<2>(indices_potentials_tuple);
        const auto contact_pot = std::get<3>(indices_potentials_tuple);

        BOOST_TEST(angle1_pot.v0()  == real_type(0.0), boost::test_tools::tolerance(tol));
        BOOST_TEST(angle1_pot.k()   == real_type(-10.0), boost::test_tools::tolerance(tol));
        BOOST_TEST(angle1_pot.n()   == std::int32_t(1));
        BOOST_TEST(angle2_pot.v0()  == real_type(1.0), boost::test_tools::tolerance(tol));
        BOOST_TEST(angle2_pot.k()   == real_type(-20.0), boost::test_tools::tolerance(tol));
        BOOST_TEST(angle2_pot.n()   == std::int32_t(2));
        BOOST_TEST(contact_pot.k()  == real_type(1.0), boost::test_tools::tolerance(tol));
        BOOST_TEST(contact_pot.v0() == real_type(2.0), boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(read_directional_contact_cosine_gaussian)
{
    mjolnir::LoggerManager::set_default_logger("test_read_directional_contact_interaction.log");
    using traits_type = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type = traits_type::real_type;
    const real_type tol = 1e-7;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
          interaction        = "DirectionalContact"
          potentials.angle1   = "Cosine"
          potentials.angle2   = "Cosine"
          potentials.contact  = "Gaussian"
          topology           = "none"
          parameters = [
          {indices = [0, 1, 2, 3], angle1 = {v0 = 0.0, k = -10.0, n = 1}, angle2 = {v0 = 0.0, k = -20.0, n = 2}, contact = {v0 = 1.0, k = 1.0, sigma = 5.0}},
          ]
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        const auto derv =
                dynamic_cast<
                    mjolnir::DirectionalContactInteraction<
                        traits_type, mjolnir::CosinePotential<real_type>,
                        mjolnir::CosinePotential<real_type>,
                        mjolnir::GaussianPotential<real_type>
                        >*
                >(base.get()); //check the expected type is contained
        BOOST_TEST(static_cast<bool>(derv));
        const auto indices_potentials_tuple = derv->potentials().at(0);
        const auto angle1_pot  = std::get<1>(indices_potentials_tuple);
        const auto angle2_pot  = std::get<2>(indices_potentials_tuple);
        const auto contact_pot = std::get<3>(indices_potentials_tuple);

        BOOST_TEST(angle1_pot.v0()  == real_type(0.0), boost::test_tools::tolerance(tol));
        BOOST_TEST(angle1_pot.k()   == real_type(-10.0), boost::test_tools::tolerance(tol));
        BOOST_TEST(angle1_pot.n()   == std::int32_t(1));
        BOOST_TEST(angle2_pot.v0()  == real_type(0.0), boost::test_tools::tolerance(tol));
        BOOST_TEST(angle2_pot.k()   == real_type(-20.0), boost::test_tools::tolerance(tol));
        BOOST_TEST(angle2_pot.n()   == std::int32_t(2));
        BOOST_TEST(contact_pot.k()  == real_type(1.0), boost::test_tools::tolerance(tol));
        BOOST_TEST(contact_pot.v0() == real_type(1.0), boost::test_tools::tolerance(tol));
        BOOST_TEST(contact_pot.sigma() == real_type(5.0), boost::test_tools::tolerance(tol));
    }
}
