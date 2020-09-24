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

    explicit RandomNumberGenerator(const std::string& internal_state)
        : nrm_(0.0, 1.0)
    {
        std::istringstream iss(internal_state);
        iss >> seed_ >> rng_;
        if(iss.fail())
        {
            throw_exception<std::runtime_error>("[error] mjolnir::"
                "RandomNumberGenerator: parse error in mt19937 ", internal_state);
        }
    }
    std::string internal_state() const
    {
        std::ostringstream oss;
        oss << seed_ << ' ' << rng_;
        return oss.str();
    }

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

    std::uint32_t seed() const noexcept {return seed_;}

  private:
    std::uint32_t seed_;
    std::mt19937        rng_;
    std::normal_distribution<real_type> nrm_;

#ifdef MJOLNIR_WITH_OPENMP
    // OpenMP implementation uses its own specialization to avoid data race.
    // So this implementation should not be instanciated with the OpenMP traits.
    static_assert(!is_openmp_simulator_traits<traits_type>::value,
                  "this is the default implementation, not for OpenMP");
#endif
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class RandomNumberGenerator<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class RandomNumberGenerator<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class RandomNumberGenerator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class RandomNumberGenerator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif /*MJOLNIR_CORE_RANDOM_NUMBER_GENERATOR*/
