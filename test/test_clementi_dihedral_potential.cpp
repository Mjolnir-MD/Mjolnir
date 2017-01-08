#define BOOST_TEST_MODULE "test_clementi_dihedral_potential"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/potential/ClementiDihedralPotential.hpp>
#include <mjolnir/core/DefaultTraits.hpp>
#include <mjolnir/util/make_unique.hpp>

typedef mjolnir::DefaultTraits traits;

constexpr static std::size_t       N    = 10000;
constexpr static traits::real_type h    = 1e-6;

BOOST_AUTO_TEST_CASE(ClementiDihedral_constructable)
{
    const traits::real_type k1 = 1.0;
    const traits::real_type k3 = 2.0;
    const traits::real_type r0 = M_PI * 2. / 3.;
    std::unique_ptr<mjolnir::LocalPotentialBase<traits>> c = 
        mjolnir::make_unique<mjolnir::ClementiDihedralPotential<traits>>(k1, k3, r0);

    BOOST_CHECK(c);
}
BOOST_AUTO_TEST_CASE(ClementiDihedral_derivative)
{
    const traits::real_type k1 = 1.0;
    const traits::real_type k3 = 2.0;
    const traits::real_type r0 = M_PI * 2. / 3.;

    auto c = mjolnir::make_unique<mjolnir::ClementiDihedralPotential<traits>>(k1, k3, r0);

    const traits::real_type x_min = 0.;
    const traits::real_type x_max = 2. * M_PI;
    const traits::real_type dx = (x_max - x_min) / N;

    traits::real_type x = x_min;
    for(std::size_t i=0; i<N; ++i)
    {
        const traits::real_type x_up = (x+h < 2. * M_PI) ? x+h : x+h - 2. * M_PI;
        const traits::real_type x_bt = (x-h > 0) ? x-h : x-h + 2. * M_PI;
        const traits::real_type pot1 = c->potential(x_up);
        const traits::real_type pot2 = c->potential(x_bt);
        const traits::real_type dpot = (pot1 - pot2) / (2 * h);
        const traits::real_type deri = c->derivative(x);

        BOOST_CHECK_CLOSE_FRACTION(dpot, deri, h);
        x += dx;
    }
}
