#define BOOST_TEST_MODULE "test_axis_aligned_plane_shape"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/AxisAlignedPlane.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <random>
#include <cmath>

BOOST_AUTO_TEST_CASE(AxisAlignedPlane_PositiveZ_geometry_unlimited)
{
    mjolnir::LoggerManager::set_default_logger("test_AxisAlignedPlane");
    using traits_type   = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type     = typename traits_type::real_type;
    using coord_type    = typename traits_type::coordinate_type;
    using boundary_type = typename traits_type::boundary_type;

    constexpr real_type tol = 1e-8;

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<double> uni(-100, 100);

    boundary_type bdry;
    mjolnir::AxisAlignedPlane<traits_type, mjolnir::PositiveZDirection> xyplane(0.0);

    for(std::size_t i=0; i<100000; ++i)
    {
        const coord_type pos(uni(mt), uni(mt), uni(mt));
        const real_type  dist = xyplane.calc_distance(pos, bdry);
        const coord_type fdir = xyplane.calc_force_direction(pos, bdry);

        BOOST_TEST(fdir[0] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(fdir[1] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(fdir[2] == 1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(dist == pos[2], boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(AxisAlignedPlane_NegativeZ_geometry_unlimited)
{
    mjolnir::LoggerManager::set_default_logger("test_AxisAlignedPlane");
    using traits_type   = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type     = typename traits_type::real_type;
    using coord_type    = typename traits_type::coordinate_type;
    using boundary_type = typename traits_type::boundary_type;

    constexpr real_type tol = 1e-8;

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<double> uni(-100, 100);

    boundary_type bdry;
    mjolnir::AxisAlignedPlane<traits_type, mjolnir::NegativeZDirection> xyplane(0.0);

    for(std::size_t i=0; i<100000; ++i)
    {
        const coord_type pos(uni(mt), uni(mt), uni(mt));
        const real_type  dist = xyplane.calc_distance(pos, bdry);
        const coord_type fdir = xyplane.calc_force_direction(pos, bdry);

        BOOST_TEST(fdir[0] ==  0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(fdir[1] ==  0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(fdir[2] == -1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(dist == -pos[2], boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(AxisAlignedPlane_geometry_periodic)
{
    mjolnir::LoggerManager::set_default_logger("test_AxisAlignedPlane");
    using traits_type   = mjolnir::SimulatorTraits<double, mjolnir::CuboidalPeriodicBoundary>;
    using real_type     = typename traits_type::real_type;
    using coord_type    = typename traits_type::coordinate_type;
    using boundary_type = typename traits_type::boundary_type;

    constexpr real_type tol = 1e-8;

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<double> uni(0, 100);

    const coord_type lower(0, 0, 0), upper(100, 100, 100);
    boundary_type bdry(lower, upper);
    mjolnir::AxisAlignedPlane<traits_type, mjolnir::PositiveZDirection> xyplane(0.0);

    for(std::size_t i=0; i<100000; ++i)
    {
        const coord_type pos(uni(mt), uni(mt), uni(mt));
        const real_type  dist = xyplane.calc_distance(pos, bdry);
        const coord_type fdir = xyplane.calc_force_direction(pos, bdry);

        const real_type d    = (pos[2] <= 50.0) ? pos[2] : pos[2] - 100.0;

        BOOST_TEST(fdir[0] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(fdir[1] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(fdir[2] == 1.0, boost::test_tools::tolerance(tol));

        if(std::abs(pos[2] - 50.0) < tol)
        {
            BOOST_TEST(std::abs(dist) == 50.0, boost::test_tools::tolerance(tol));
        }
        else
        {
            BOOST_TEST(dist == d, boost::test_tools::tolerance(tol));
        }
    }
}

template<typename realT>
struct dummy_potential
{
    realT cutoff;
    std::size_t N;

    std::vector<std::size_t> participants() const
    {
        std::vector<std::size_t> retval(N);
        std::iota(retval.begin(), retval.end(), 0);
        return retval;
    }

    realT max_cutoff_length() const noexcept {return cutoff;}

};

BOOST_AUTO_TEST_CASE(AxisAlignedPlane_neighbors_unlimited)
{
    mjolnir::LoggerManager::set_default_logger("test_AxisAlignedPlane");
    using traits_type   = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using coord_type    = typename traits_type::coordinate_type;
    using boundary_type = typename traits_type::boundary_type;
    using system_type   = mjolnir::System<traits_type>;

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<double> uni(-10, 10);

    system_type sys(1000, boundary_type{});
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        auto p = sys[i];
        p.mass     = 0.0;
        p.position = coord_type(uni(mt), uni(mt), uni(mt));
        p.velocity = coord_type(0, 0, 0);
        p.force    = coord_type(0, 0, 0);
    }

    dummy_potential<double> dummy; dummy.cutoff = 1.0; dummy.N = 1000;
    mjolnir::AxisAlignedPlane<traits_type, mjolnir::PositiveZDirection>
        xyplane(0.0, 0.0); // no mergin here
    xyplane.initialize(sys, dummy);

    const auto neighbors = xyplane.neighbors();

    for(std::size_t i : neighbors)
    {
        if(std::abs(sys.at(i).position[2]) > 1.0)
        {
            std::cerr << i << ", " << sys.at(i).position[2] << std::endl;
        }
        BOOST_TEST(std::abs(sys.at(i).position[2]) <= 1.0);
    }

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        if(std::abs(sys.at(i).position[2]) > 1.0){continue;}
        const bool found =
            std::find(neighbors.begin(), neighbors.end(), i) != neighbors.end();
        BOOST_TEST(found);
    }
}
