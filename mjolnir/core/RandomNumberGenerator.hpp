#ifndef MJOLNIR_CORE_RANDOM_NUMBER_GENERATOR_HPP
#define MJOLNIR_CORE_RANDOM_NUMBER_GENERATOR_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <random>
#include <cstdint>

namespace mjolnir
{

template<typename traitsT>
class RandomNumberGenerator
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;

  public:
    explicit RandomNumberGenerator(const std::uint32_t seed)
        : seed_(seed), rng_(seed), nrm_(0.0, 1.0)
    {}
    ~RandomNumberGenerator() = default;

    real_type uniform_real01()
    {
        return std::generate_canonical<
            real_type, std::numeric_limits<real_type>::digits
            >(this->rng_);
    }
    real_type uniform_real(const real_type min, const real_type max)
    {
        return this->uniform_real01() * (max - min) + min;
    }

    real_type gaussian()
    {
        return this->nrm_(this->rng_);
    }
    real_type gaussian(const real_type mean, const real_type stddev)
    {
        return this->nrm_(this->rng_) * stddev + mean;
    }

  private:
    const std::uint32_t seed_;
    std::mt19937        rng_;
    std::normal_distribution<real_type> nrm_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class RandomNumberGenerator<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class RandomNumberGenerator<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class RandomNumberGenerator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class RandomNumberGenerator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif /*MJOLNIR_CORE_RANDOM_NUMBER_GENERATOR*/
