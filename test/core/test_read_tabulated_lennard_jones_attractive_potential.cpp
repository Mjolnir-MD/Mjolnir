#define BOOST_TEST_MODULE "test_read_tabulated_lennard_jones_attractive_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/input/read_global_potential.hpp>
#include <tuple>

using test_types = std::tuple<double, float>;

constexpr inline float  tolerance_value(float)  noexcept {return 1e-4;}
constexpr inline double tolerance_value(double) noexcept {return 1e-8;}

template<typename Real>
decltype(boost::test_tools::tolerance(std::declval<Real>()))
tolerance() {return boost::test_tools::tolerance(tolerance_value(Real()));}

BOOST_AUTO_TEST_CASE_TEMPLATE(read_tabulated_lennard_jones_attractive_noenv, T, test_types)
{
    mjolnir::LoggerManager::set_default_logger("test_read_tabulated_lennard_jones_attractive.log");

    using real_type = T;
    using traits_type = mjolnir::SimulatorTraits<real_type, mjolnir::UnlimitedBoundary>;
    using potential_type = mjolnir::LennardJonesAttractivePotential<real_type>;
    {
        using namespace toml::literals;
        const auto v = u8R"(
            interaction = "Pair"
            potential   = "LennardJonesAttractive"
            spatial_partition.type  = "Naive"
            ignore.molecule         = "Nothing"
            ignore.particles_within.bond    = 3
            ignore.particles_within.contact = 1
            table.A.A = {sigma = 1.0, epsilon = 0.5}
            table.A.B = {sigma = 2.0, epsilon = 1.5}
            table.B.B = {sigma = 3.0, epsilon = 2.5}
            parameters = [
                {index =   0, name = "A"},
                {index =   1, name = "B"},
                {index =   2, name = "A"},
                {index =   3, name = "B"},
                {index =   5, name = "A"},
                {index =   7, name = "B"},
                {index = 100, name = "A"},
            ]
        )"_toml;

        const auto pot_para = mjolnir::read_lennard_jones_attractive_potential<traits_type>(v);
        const auto& para = dynamic_cast<
            mjolnir::CombinationTable<traits_type, potential_type> const&
            >(pot_para.second.cref());

        const auto ignore_within = para.exclusion_list().ignore_topology();
        const std::map<std::string, std::size_t> within(
                ignore_within.begin(), ignore_within.end());

        BOOST_TEST(within.size() == 2u);
        BOOST_TEST(within.at("bond")    == 3ul);
        BOOST_TEST(within.at("contact") == 1ul);

        BOOST_TEST(!para.exclusion_list().is_ignored_molecule(0, 0));
        BOOST_TEST(!para.exclusion_list().is_ignored_molecule(0, 1));
        BOOST_TEST(!para.exclusion_list().is_ignored_molecule(1, 1));

        BOOST_TEST(para.participants().size() ==   7u);
        BOOST_TEST(para.participants().at(0)  ==   0u);
        BOOST_TEST(para.participants().at(1)  ==   1u);
        BOOST_TEST(para.participants().at(2)  ==   2u);
        BOOST_TEST(para.participants().at(3)  ==   3u);
        BOOST_TEST(para.participants().at(4)  ==   5u);
        BOOST_TEST(para.participants().at(5)  ==   7u);
        BOOST_TEST(para.participants().at(6)  == 100u);

        BOOST_TEST(para.parameters().at(  0)  == "A", tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  1)  == "B", tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  2)  == "A", tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  3)  == "B", tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  5)  == "A", tolerance<real_type>());
        BOOST_TEST(para.parameters().at(  7)  == "B", tolerance<real_type>());
        BOOST_TEST(para.parameters().at(100)  == "A", tolerance<real_type>());

        const auto para_AA = para.prepare_params(0, 2);
        const auto para_AB = para.prepare_params(0, 1);
        const auto para_BA = para.prepare_params(1, 2);
        const auto para_BB = para.prepare_params(1, 3);

        BOOST_TEST(para_AA.sigma   == 1.0, tolerance<real_type>());
        BOOST_TEST(para_AA.epsilon == 0.5, tolerance<real_type>());
        BOOST_TEST(para_AB.sigma   == 2.0, tolerance<real_type>());
        BOOST_TEST(para_AB.epsilon == 1.5, tolerance<real_type>());
        BOOST_TEST(para_BA.sigma   == 2.0, tolerance<real_type>());
        BOOST_TEST(para_BA.epsilon == 1.5, tolerance<real_type>());
        BOOST_TEST(para_BB.sigma   == 3.0, tolerance<real_type>());
        BOOST_TEST(para_BB.epsilon == 2.5, tolerance<real_type>());
    }
}
