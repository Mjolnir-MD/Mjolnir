#define BOOST_TEST_MODULE "test_excluded_volume_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/forcefield/external/ExcludedVolumeWallPotential.hpp>

BOOST_AUTO_TEST_CASE(ExcludedVolumeWallPotential_double)
{
    using real_type = double;
    constexpr static std::size_t N = 10000;
    constexpr static real_type h   = 1e-6;
    constexpr static real_type tol = 1e-5;

    const real_type epsilon = 1.0;
    const real_type radius  = 1.0;
    const std::vector<std::pair<std::size_t, real_type>> radii{{0, radius}};

    mjolnir::ExcludedVolumeWallPotential<real_type> exvw(epsilon, 2.0, radii);

    const real_type cutoff_length = exvw.max_cutoff_length();
    const real_type z_max         = cutoff_length;
    const real_type z_min         = radius * 0.5;
    const real_type dz            = (z_max - z_min) / N;

    real_type z = z_min;
    for(std::size_t i = 0; i < N; ++i)
    {
        const real_type pot1 = exvw.potential(0, z + h);
        const real_type pot2 = exvw.potential(0, z - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = exvw.derivative(0, z);

        if(std::abs(deri) > tol)
        {
            BOOST_TEST(dpot == deri, boost::test_tools::tolerance(tol));
        }
        else
        {
            BOOST_TEST(deri == 0.0, boost::test_tools::tolerance(tol));
        }
        z += dz;
    }
}
