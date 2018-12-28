#define BOOST_TEST_MODULE "test_read_forcefield"

#include <boost/test/included/unit_test.hpp>
#include <mjolnir/input/read_forcefield.hpp>

BOOST_AUTO_TEST_CASE(read_empty_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    constexpr real_type tol = 1e-8;
    {
        const toml::table v = toml::table{
            {"local",       toml::array{}},
            {"global",      toml::array{}},
            {"external",    toml::array{}}
        };

        const auto ff = mjolnir::read_forcefield_from_table<traits_type>(v, "./");

        BOOST_TEST(ff.local().size()    == 0);
        BOOST_TEST(ff.global().size()   == 0);
        BOOST_TEST(ff.external().size() == 0);
    }
}

BOOST_AUTO_TEST_CASE(read_several_forcefield)
{
    mjolnir::LoggerManager::set_default_logger("test_read_forcefield.log");

    using real_type = double;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    constexpr real_type tol = 1e-8;
    {
        const toml::table v = toml::table{
            {"local",       toml::array{
                toml::table{
                    {"interaction", toml::value("BondAngle")},
                    {"potential",   toml::value("Harmonic")},
                    {"topology",    toml::value("none")},
                    {"parameters",  toml::value(toml::array(/*empty*/))}
                }, toml::table{
                    {"interaction", toml::value("BondLength")},
                    {"potential",   toml::value("Harmonic")},
                    {"topology",    toml::value("bond")},
                    {"parameters",  toml::value(toml::array(/*empty*/))}
                }
            }},
            {"global",      toml::array{
                toml::table{
                    {"interaction",       toml::value("Pair")},
                    {"potential",         toml::value("ExcludedVolume")},
                    {"spatial_partition", toml::value(toml::table{
                                {"type", toml::value("Naive")}
                    })},
                    {"epsilon",           toml::value(3.14)},
                    {"ignore",            toml::value(toml::table{
                        {"molecule",         toml::value("Nothing")},
                        {"particles_within", toml::table{{"bond", 3}, {"contact", 1}}},
                    })},
                    {"parameters",        toml::value(toml::array())}
                }, toml::table{
                    {"interaction",       toml::value("Pair")},
                    {"potential",         toml::value("LennardJones")},
                    {"spatial_partition", toml::value(toml::table{
                                {"type", toml::value("Naive")}
                    })},
                    {"ignore",            toml::value(toml::table{
                        {"molecule",         toml::value("Nothing")},
                        {"particles_within", toml::table{{"bond", 3}, {"contact", 1}}},
                    })},
                    {"parameters",        toml::value(toml::array{})}
                }
            }},
            {"external",    toml::array{
                toml::table{
                    {"interaction",       toml::value("Distance")},
                    {"potential",         toml::value("ExcludedVolumeWall")},
                    {"shape", toml::value(toml::table{
                            {"name",     toml::value("AxisAlignedPlane")},
                            {"axis",     toml::value("+X")},
                            {"position", toml::value(1.0)},
                            {"margin",   toml::value(0.5)}
                    })},
                    {"epsilon",           toml::value(3.14)},
                    {"parameters",        toml::value(toml::array{})}
                }, toml::table{
                    {"interaction",       toml::value("Distance")},
                    {"potential",         toml::value("LennardJonesWall")},
                    {"shape", toml::value(toml::table{
                            {"name",     toml::value("AxisAlignedPlane")},
                            {"axis",     toml::value("+X")},
                            {"position", toml::value(1.0)},
                            {"margin",   toml::value(0.5)}
                    })},
                    {"parameters",        toml::value(toml::array{})}
                }
            }}
        };

        const auto ff = mjolnir::read_forcefield_from_table<traits_type>(v, "./");
        BOOST_TEST(ff.local().size()    == 2);
        BOOST_TEST(ff.global().size()   == 2);
        BOOST_TEST(ff.external().size() == 2);

        // contents are tested in the individual test codes
    }
}
