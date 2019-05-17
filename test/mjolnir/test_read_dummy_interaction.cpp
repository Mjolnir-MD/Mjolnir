#define BOOST_TEST_MODULE "test_read_dummy_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_local_interaction.hpp>
#include <test/util/traits.hpp>

BOOST_AUTO_TEST_CASE(read_dummy_interaction)
{
    mjolnir::LoggerManager::set_default_logger("test_read_dummy_interaction.log");

    using real_type = double;
    using traits_type = mjolnir::test::traits<real_type>;
    {
        using namespace toml::literals;
        const toml::value v = u8R"(
            interaction = "Dummy"
            topology    = "bond"
            parameters  = [
                {indices = [0, 1]},
                {indices = [1, 2]},
            ] # empty
        )"_toml;

        const auto base = mjolnir::read_local_interaction<traits_type>(v);
        BOOST_TEST(static_cast<bool>(base));

        // check the expected type is contained
        const auto derv = dynamic_cast<mjolnir::DummyInteraction<
            traits_type>*>(base.get());
        BOOST_TEST(static_cast<bool>(derv));

        // check the correct topology is contained

        mjolnir::Topology topol(3);
        derv->write_topology(topol);

        const auto adj_0 = topol.list_adjacent_within(0, 1, "bond"); // [1,0]
        const auto adj_1 = topol.list_adjacent_within(1, 1, "bond"); // [0,1,2]
        const auto adj_2 = topol.list_adjacent_within(2, 1, "bond"); // [1,2]

        BOOST_TEST(adj_0.size() == 2u);
        BOOST_TEST(adj_0.at(0)  == 0u);
        BOOST_TEST(adj_0.at(1)  == 1u);

        BOOST_TEST(adj_1.size() == 3u);
        BOOST_TEST(adj_1.at(0)  == 0u);
        BOOST_TEST(adj_1.at(1)  == 1u);
        BOOST_TEST(adj_1.at(2)  == 2u);

        BOOST_TEST(adj_2.size() == 2u);
        BOOST_TEST(adj_2.at(0)  == 1u);
        BOOST_TEST(adj_2.at(1)  == 2u);


    }
}
