#define BOOST_TEST_MODULE "test_implicit_membrane_potential"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/forcefield/external/ImplicitMembranePotential.hpp>

BOOST_AUTO_TEST_CASE(ImplicitMembranePotential_double)
{
    using real_type = double;
    constexpr static std::size_t N = 10000;
    constexpr static real_type h   = 1e-6;
    constexpr static real_type tol = 1e-5;

    const real_type thickness = 10.0;
    const real_type magnitude = 1.0;
    const real_type bend = 1.5;
    const std::vector<std::pair<std::size_t, real_type>> hydrophobicities{
        {0, 1.}, {1, 0.}
    };

    mjolnir::ImplicitMembranePotential<real_type>
        im(thickness, magnitude, bend, 4.0, hydrophobicities);

    const real_type cutoff_length = im.max_cutoff_length();
    const real_type z_min = -1 * cutoff_length;
    const real_type z_max = cutoff_length;
    const real_type dz = (z_max - z_min) / N;

    real_type z = z_min;
    for(std::size_t i = 0; i < N; ++i)
    {
        const real_type pot1 = im.potential(0, z + h);
        const real_type pot2 = im.potential(0, z - h);
        const real_type dpot = (pot1 - pot2) / (2 * h);
        const real_type deri = im.derivative(0, z);

        if(std::abs(z) > h)
        {
            if(std::abs(deri) > tol)
            {
                BOOST_TEST(dpot == deri, boost::test_tools::tolerance(tol));
            }
            else
            {
                BOOST_TEST(deri == 0.0, boost::test_tools::tolerance(tol));
            }
        }

        const real_type pot0 = im.potential(1, z);
        const real_type deri0 = im.derivative(1, z);

        BOOST_TEST(pot0  == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(deri0 == 0.0, boost::test_tools::tolerance(tol));

        z += dz;
    }
}
