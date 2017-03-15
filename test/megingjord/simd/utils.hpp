#ifndef MJOLNIR_TEST_UTILS
#define MJOLNIR_TEST_UTILS

#ifndef BOOST_TEST_MODULE
#error "this header file is for test. DO NOT USE in normal code"
#endif

#include <random>

constexpr static unsigned int seed = 32479327;
constexpr static std::size_t N = 10000;
constexpr static double tolerance = 1e-12;

constexpr static std::size_t tight = 8;
constexpr static std::size_t loose = 9;

template<typename T, std::size_t N>
struct array_maker
{
    template<typename gen, typename distro>
    static std::array<T, N> invoke(gen& g, distro& dist)
    {
        std::array<T, N> retval;
        for(std::size_t i=0; i<N; ++i)
            retval[i] = dist(g);
        return retval;
    }
};
#endif// MJOLNIR_TEST_UTILS
