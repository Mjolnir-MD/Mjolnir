#define BOOST_TEST_MODULE "test_axis_aligned_plane_shape"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <test/mjolnir/traits.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/AxisAlignedPlane.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <random>
#include <cmath>

BOOST_AUTO_TEST_CASE(AxisAlignedPlane_geometry_unlimited)
{
    using traits = mjolnir::test::traits<double>;
    using real_type     = typename traits::real_type;
    using coord_type    = typename traits::coordinate_type;
    using boundary_type = typename traits::boundary_type;

    static constexpr traits::real_type tolerance = 1e-8;

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<double> uni(-100, 100);

    boundary_type bdry;
    mjolnir::AxisAlignedPlane<traits, mjolnir::PositiveZDirection> xyplane(0.0);

    for(std::size_t i=0; i<100000; ++i)
    {
        const coord_type pos(uni(mt), uni(mt), uni(mt));
        const real_type  dist = xyplane.calc_distance(pos, bdry);
        const coord_type fdir = xyplane.calc_force_direction(pos, bdry);
        const real_type  sign = std::copysign(1.0, pos[2]);

        BOOST_CHECK_SMALL(fdir[0], tolerance);
        BOOST_CHECK_SMALL(fdir[1], tolerance);
        BOOST_CHECK_EQUAL(fdir[2], sign);
        BOOST_CHECK_CLOSE_FRACTION(dist, pos[2], tolerance);
    }
}

BOOST_AUTO_TEST_CASE(AxisAlignedPlane_geometry_periodic)
{
    using traits = mjolnir::test::traits<double, mjolnir::CubicPeriodicBoundary>;
    using real_type     = typename traits::real_type;
    using coord_type    = typename traits::coordinate_type;
    using boundary_type = typename traits::boundary_type;

    static constexpr traits::real_type tolerance = 1e-8;

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<double> uni(0, 100);

    const coord_type lower(0, 0, 0), upper(100, 100, 100);
    boundary_type bdry(lower, upper);
    mjolnir::AxisAlignedPlane<traits, mjolnir::PositiveZDirection> xyplane(0.0);

    for(std::size_t i=0; i<100000; ++i)
    {
        const coord_type pos(uni(mt), uni(mt), uni(mt));
        const real_type  dist = xyplane.calc_distance(pos, bdry);
        const coord_type fdir = xyplane.calc_force_direction(pos, bdry);

        const real_type sign = (pos[2] <= 50.0) ? 1.0    : -1.0;
        const real_type d    = (pos[2] <= 50.0) ? pos[2] : pos[2] - 100.0;

        BOOST_CHECK_SMALL(fdir[0], tolerance);
        BOOST_CHECK_SMALL(fdir[1], tolerance);
        BOOST_CHECK_EQUAL(fdir[2], sign);

        if(std::abs(pos[2] - 50.0) < tolerance)
        {
            BOOST_CHECK_CLOSE_FRACTION(std::abs(dist), 50.0, tolerance);
        }
        else
        {
            BOOST_CHECK_CLOSE_FRACTION(dist, d, tolerance);
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
    using traits = mjolnir::test::traits<double>;
    using real_type     = typename traits::real_type;
    using coord_type    = typename traits::coordinate_type;
    using boundary_type = typename traits::boundary_type;
    using system_type   = mjolnir::System<traits>;
    using particle_type = typename system_type::particle_type;

    static constexpr traits::real_type tolerance = 1e-8;

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<double> uni(-10, 10);

    std::vector<particle_type> pcon(1000);
    for(auto& p : pcon)
    {
        p.mass     = 0.0;
        p.position = coord_type(uni(mt), uni(mt), uni(mt));
        p.velocity = coord_type(0, 0, 0);
        p.force    = coord_type(0, 0, 0);
    }
    boundary_type bdry;
    system_type sys(std::move(pcon), std::move(bdry));

    dummy_potential<double> dummy; dummy.cutoff = 1.0; dummy.N = 1000;
    mjolnir::AxisAlignedPlane<traits, mjolnir::PositiveZDirection> xyplane(0.0, 0.0); // no mergin here
    xyplane.initialize(sys, dummy);

    const auto neighbors = xyplane.neighbors();

    for(std::size_t i : neighbors)
    {
        if(std::abs(sys.at(i).position[2]) > 1.0)
        {
            std::cerr << i << ", " << sys.at(i).position[2] << std::endl;
        }
        BOOST_CHECK(std::abs(sys.at(i).position[2]) <= 1.0);
    }

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        if(std::abs(sys.at(i).position[2]) > 1.0){continue;}
        BOOST_CHECK(std::find(neighbors.begin(), neighbors.end(), i)
                    != neighbors.end());
    }
}
