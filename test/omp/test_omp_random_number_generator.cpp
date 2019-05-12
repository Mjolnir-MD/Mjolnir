#define BOOST_TEST_MODULE "test_omp_random_number_generator"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/omp/RandomNumberGenerator.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <random>
#include <algorithm>


// It tests that ...
// - 
//
//
BOOST_AUTO_TEST_CASE(test_omp_random_number_generator)
{
    BOOST_TEST_WARN(omp_get_max_threads() > 2);
    BOOST_TEST_MESSAGE("maximum number of threads = " << omp_get_max_threads());

    mjolnir::LoggerManager::set_default_logger("test_omp_random_number_generator.log");

    const std::uint32_t seed=123456789;

    mjolnir::RandomNumberGenerator<
        mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>
        > rng_1(seed);
    mjolnir::RandomNumberGenerator<
        mjolnir::OpenMPSimulatorTraits<double, mjolnir::UnlimitedBoundary>
        > rng_2(seed);

#pragma omp parallel for
    for(std::size_t i=0; i<10000; ++i)
    {
        // this part is called in parallel
        // {
        const auto real01_1 = rng_1.uniform_real01();
        const auto real01_2 = rng_2.uniform_real01();
        // }

#pragma omp critical
        {
            BOOST_TEST(0.0 <= real01_1);
            BOOST_TEST(real01_1 < 1.0);
            BOOST_TEST(0.0 <= real01_2);
            BOOST_TEST(real01_2 < 1.0);
            BOOST_TEST(real01_1 == real01_2);
        }
    }
}
