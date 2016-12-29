#ifndef MJOLNIR_CORE_RANDOM_NUMBER_GENERATOR
#define MJOLNIR_CORE_RANDOM_NUMBER_GENERATOR
#include <random>

namespace mjolnir
{

template<typename traitsT>
class RandomNumberGenerator
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:
    RandomNumberGenerator(const std::uint_fast32_t seed)
        : seed_(seed), rng_(seed)
    {}
    ~RandomNumberGenerator() = default;

    real_type uniform_real(const real_type min, const real_type max);
    real_type gaussian(const real_type mean, const real_type sigma);

    coordinate_type
    maxwell_boltzmann(const real_type mass, const real_type temperature,
                      const real_type kB);

    coordinate_type
    underdamped_langevin(const real_type mass, const real_type gamma,
        const time_type dt, const real_type temperature, const real_type kB);

  private:
    const std::uint_fast32_t seed_;
    std::mt19937 rng_;
};

// TODO: measure
// mt19937 v.s. mt19937_64
// (or other random number generator like xorshift)
//
// (construct distribution at each time)
// v.s.
// (pre-constructed distro and some calculation)
//
// and
//
// standard C++ library
// v.s.
// GNU Scientific library [[likely]]
 

template<typename traitsT>
inline typename RandomNumberGenerator<traitsT>::real_type
RandomNumberGenerator<traitsT>::uniform_real(
    const real_type min, const real_type max)
{
    return (std::uniform_real_distribution<real_type>(min, max))(rng_);
}

template<typename traitsT>
inline typename RandomNumberGenerator<traitsT>::real_type
RandomNumberGenerator<traitsT>::gaussian(
    const real_type mean, const real_type sigma)
{
    return (std::normal_distribution<real_type>(mean, sigma))(rng_);
}

template<typename traitsT>
inline typename RandomNumberGenerator<traitsT>::coordinate_type
RandomNumberGenerator<traitsT>::maxwell_boltzmann(
    const real_type m, const real_type T, const real_type kB)
{
    std::normal_distribution<real_type> distro(0., std::sqrt(kB * T / m));
    return coordinate_type(distro(rng_), distro(rng_), distro(rng_));
}

template<typename traitsT>
inline typename RandomNumberGenerator<traitsT>::coordinate_type
RandomNumberGenerator<traitsT>::underdamped_langevin(const real_type m,
        const real_type g, const time_type dt, const real_type T, const real_type kB)
{
    std::normal_distribution<real_type>
        distro(0., std::sqrt(2 * kB * T * g / (m * dt)));
    return coordinate_type(distro(rng_), distro(rng_), distro(rng_));
}

} // mjolnir
#endif /*MJOLNIR_CORE_RANDOM_NUMBER_GENERATOR*/
