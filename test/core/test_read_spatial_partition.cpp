#define BOOST_TEST_MODULE "test_read_spatial_partition"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/forcefield/global/DebyeHuckelPotential.hpp>
#include <mjolnir/input/read_spatial_partition.hpp>

BOOST_AUTO_TEST_CASE(read_spatial_partition_naive)
{
    mjolnir::LoggerManager::set_default_logger("test_read_spatial_partition.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::DebyeHuckelPotential<traits_type>;

    {
        using namespace toml::literals;
        const toml::table v = toml::get<toml::table>(u8R"(
            spatial_partition.type = "Naive"
        )"_toml);

        const auto sp =
            mjolnir::read_spatial_partition<traits_type, potential_type>(v);

        const auto sp_ptr = std::addressof(sp.base());
        BOOST_TEST(sp_ptr);

        const auto naive_ptr = dynamic_cast<
            mjolnir::NaivePairCalculation<traits_type, potential_type> const*
            >(sp_ptr);
        BOOST_TEST(naive_ptr);
    }
}

BOOST_AUTO_TEST_CASE(read_spatial_partition_verlet)
{
    mjolnir::LoggerManager::set_default_logger("test_read_spatial_partition.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::DebyeHuckelPotential<traits_type>;

    {
        using namespace toml::literals;
        const toml::table v = toml::get<toml::table>(u8R"(
            spatial_partition.type = "VerletList"
            spatial_partition.margin = 0.5
        )"_toml);

        const auto sp =
            mjolnir::read_spatial_partition<traits_type, potential_type>(v);
        BOOST_TEST(sp.margin() == 0.5);

        const auto sp_ptr = std::addressof(sp.base());
        BOOST_TEST(sp_ptr);

        const auto verlet_ptr = dynamic_cast<
            mjolnir::VerletList<traits_type, potential_type> const*
            >(sp_ptr);
        BOOST_TEST(verlet_ptr);
    }
}

BOOST_AUTO_TEST_CASE(read_spatial_partition_unlimited_cell_list)
{
    mjolnir::LoggerManager::set_default_logger("test_read_spatial_partition.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::DebyeHuckelPotential<traits_type>;

    {
        using namespace toml::literals;
        const toml::table v = toml::get<toml::table>(u8R"(
            spatial_partition.type = "CellList"
            spatial_partition.margin = 0.5
        )"_toml);

        const auto sp =
            mjolnir::read_spatial_partition<traits_type, potential_type>(v);
        BOOST_TEST(sp.margin() == 0.5);

        const auto sp_ptr = std::addressof(sp.base());
        BOOST_TEST(sp_ptr);

        const auto cell_ptr = dynamic_cast<
            mjolnir::UnlimitedGridCellList<traits_type, potential_type> const*
            >(sp_ptr);
        BOOST_TEST(cell_ptr);
    }
}

BOOST_AUTO_TEST_CASE(read_spatial_partition_periodic_cell_list)
{
    mjolnir::LoggerManager::set_default_logger("test_read_spatial_partition.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::CuboidalPeriodicBoundary>;
    using potential_type = mjolnir::DebyeHuckelPotential<traits_type>;

    {
        using namespace toml::literals;
        const toml::table v = toml::get<toml::table>(u8R"(
            spatial_partition.type = "CellList"
            spatial_partition.margin = 0.5
        )"_toml);

        const auto sp =
            mjolnir::read_spatial_partition<traits_type, potential_type>(v);
        BOOST_TEST(sp.margin() == 0.5);

        const auto sp_ptr = std::addressof(sp.base());
        BOOST_TEST(sp_ptr);

        const auto cell_ptr = dynamic_cast<
            mjolnir::PeriodicGridCellList<traits_type, potential_type> const*
            >(sp_ptr);
        BOOST_TEST(cell_ptr);
    }
}
