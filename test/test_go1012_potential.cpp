#define BOOST_TEST_MODULE "test_go1012_potential"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/Go1012ContactPotential.hpp>
#include <mjolnir/core/DefaultTraits.hpp>
#include <mjolnir/util/make_unique.hpp>

typedef mjolnir::DefaultTraits traits;

constexpr static std::size_t       N = 10000;
constexpr static traits::real_type h = 1e-6;

BOOST_AUTO_TEST_CASE(Go1012_constructable)
{
    const traits::real_type e  = 1.0;
    const traits::real_type r0 = 5.0;
    std::unique_ptr<mjolnir::LocalPotentialBase<traits>> Go = 
        mjolnir::make_unique<mjolnir::Go1012ContactPotential<traits>>(e, r0);

    BOOST_CHECK(Go);
}
BOOST_AUTO_TEST_CASE(Harmonic_derivative)
{
    const traits::real_type e  = 1.0;
    const traits::real_type r0 = 5.0;

    auto Go =  mjolnir::make_unique<mjolnir::Go1012ContactPotential<traits>>(e, r0);

    const traits::real_type x_min = 0.8 * r0;
    const traits::real_type x_max = 5.0 * r0;
    const traits::real_type dx = (x_max - x_min) / N;

    traits::real_type x = x_min;
    for(std::size_t i=0; i<N; ++i)
    {
        const traits::real_type pot1 = Go->potential(x + h);
        const traits::real_type pot2 = Go->potential(x - h);
        const traits::real_type dpot = (pot1 - pot2) / (2 * h);
        const traits::real_type deri = Go->derivative(x);

        BOOST_CHECK_CLOSE_FRACTION(dpot, deri, h);
        x += dx;
    }
}
