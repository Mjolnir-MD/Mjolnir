#ifndef MJOLNIR_CORE_RANDOM_NUMBER_GENERATOR
#define MJOLNIR_CORE_RANDOM_NUMBER_GENERATOR
#include <random>
#include <cstdint>

namespace mjolnir
{

template<typename traitsT>
class RandomNumberGenerator
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:
    RandomNumberGenerator(const std::uint32_t seed)
        : seed_(seed), rng_(seed)
    {}
    ~RandomNumberGenerator() = default;

    real_type uniform_real(const real_type min, const real_type max)
    {
        return this->uni_(this->rng_) * (max - min) + min;
    }
    real_type gaussian(const real_type mean, const real_type stddev)
    {
        return this->nrm_(this->rng_) * stddev + mean;
    }

  private:
    const std::uint32_t seed_;
    std::mt19937        rng_;
    std::uniform_real_distribution<real_type> uni_;
    std::normal_distribution<real_type>       nrm_;
};

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

} // mjolnir
#endif /*MJOLNIR_CORE_RANDOM_NUMBER_GENERATOR*/
