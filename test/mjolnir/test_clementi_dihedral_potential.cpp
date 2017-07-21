#define BOOST_TEST_MODULE "test_clementi_dihedral_potential"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/potential/local/ClementiDihedralPotential.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/constants.hpp>
#include <mjolnir/util/make_unique.hpp>



BOOST_AUTO_TEST_CASE(ClementiDihedral_double)
{
    typedef mjolnir::SimulatorTraitsBase<double, mjolnir::UnlimitedBoundary> traits;
    constexpr static std::size_t       N   = 1000;
    constexpr static traits::real_type h   = 1e-6;
    constexpr static traits::real_type pi  = mjolnir::constants<traits::real_type>::pi;

    const traits::real_type k1 = 1.0;
    const traits::real_type k3 = 2.0;
    const traits::real_type r0 = pi * 2. / 3.;

    mjolnir::ClementiDihedralPotential<traits> c(k1, k3, r0);

    const traits::real_type x_min = 0.;
    const traits::real_type x_max = 2. * pi;
    const traits::real_type dx = (x_max - x_min) / N;

    traits::real_type x = x_min;
    for(std::size_t i=0; i<N; ++i)
    {
        const traits::real_type x_up = (x+h < 2. * pi) ? x+h : x+h - 2. * pi;
        const traits::real_type x_bt = (x-h > 0) ? x-h : x-h + 2. * pi;
        const traits::real_type pot1 = c.potential(x_up);
        const traits::real_type pot2 = c.potential(x_bt);
        const traits::real_type dpot = (pot1 - pot2) / (2 * h);
        const traits::real_type deri = c.derivative(x);

        if(std::abs(deri) > h)
            BOOST_CHECK_CLOSE_FRACTION(dpot, deri, h);
        else
            BOOST_CHECK_SMALL(deri, h);

        x += dx;
    }
}

BOOST_AUTO_TEST_CASE(ClementiDihedral_float)
{
    typedef mjolnir::SimulatorTraitsBase<float, mjolnir::UnlimitedBoundary> traits;
    constexpr static std::size_t       N    = 100;
    constexpr static traits::real_type h    = 1e-3;
    constexpr static traits::real_type pi = mjolnir::constants<traits::real_type>::pi;

    const traits::real_type k1 = 1.0;
    const traits::real_type k3 = 2.0;
    const traits::real_type r0 = pi * 2. / 3.;

    mjolnir::ClementiDihedralPotential<traits> c(k1, k3, r0);

    const traits::real_type x_min = 0.;
    const traits::real_type x_max = 2. * pi;
    const traits::real_type dx = (x_max - x_min) / N;

    traits::real_type x = x_min;
    for(std::size_t i=0; i<N; ++i)
    {
        const traits::real_type x_up = (x+h < 2. * pi) ? x+h : x+h - 2. * pi;
        const traits::real_type x_bt = (x-h > 0) ? x-h : x-h + 2. * pi;
        const traits::real_type pot1 = c.potential(x_up);
        const traits::real_type pot2 = c.potential(x_bt);
        const traits::real_type dpot = (pot1 - pot2) / (2 * h);
        const traits::real_type deri = c.derivative(x);

        if(std::abs(deri) > h)
            BOOST_CHECK_CLOSE_FRACTION(dpot, deri, h);
        else
            BOOST_CHECK_SMALL(deri, h);

        x += dx;
    }
}

