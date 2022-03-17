#ifndef MJOLNIR_POTENTIAL_GLOBAL_LENNARD_JONES_ATTRACTIVE_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_LENNARD_JONES_ATTRACTIVE_POTENTIAL_HPP
#include <algorithm>
#include <utility>
#include <cmath>

namespace mjolnir
{

template<typename T> class System;

template<typename realT>
class LennardJonesAttractivePotential
{
  public:
    using real_type = realT;

    struct parameter_type
    {
        real_type sigma;
        real_type epsilon;
    };

    static constexpr real_type default_cutoff() noexcept
    {
        return real_type(2.5);
    }

  public:

    explicit LennardJonesAttractivePotential(const real_type cutoff_ratio) noexcept
        : cutoff_ratio_(cutoff_ratio), coef_at_cutoff_(
                std::pow(real_type(1) / cutoff_ratio, 12) -
                std::pow(real_type(1) / cutoff_ratio,  6))
    {}
    ~LennardJonesAttractivePotential() = default;

    real_type potential(const real_type r, const parameter_type& params) const noexcept
    {
        constexpr real_type sixth_root_of_two(1.12246204831);

        if(r < params.sigma * sixth_root_of_two)
        {
            return -params.epsilon;
        }
        else if(params.sigma * cutoff_ratio_ < r)
        {
            return 0;
        }
        const real_type sr1 = params.sigma / r;
        const real_type sr3 = sr1 * sr1 * sr1;
        const real_type sr6 = sr3 * sr3;
        return 4 * params.epsilon * (sr6 * (sr6 - 1) - coef_at_cutoff_);
    }
    real_type derivative(const real_type r, const parameter_type& params) const noexcept
    {
        constexpr real_type sixth_root_of_two(1.12246204831);

        if(r < params.sigma * sixth_root_of_two ||
               params.sigma * cutoff_ratio_ < r)
        {
            return 0;
        }

        const real_type rinv = 1 / r;
        const real_type sr1 = params.sigma * rinv;
        const real_type sr3 = sr1 * sr1 * sr1;
        const real_type sr6 = sr3 * sr3;
        return 24 * params.epsilon * (sr6 - 2 * sr6 * sr6) * rinv;
    }

    template<typename T>
    void initialize(const System<T>&) noexcept {return;}

    template<typename T>
    void update(const System<T>&) noexcept {return;}

    static const char* name() noexcept {return "LennardJonesAttractive";}

    // It takes per-particle parameters and return the maximum cutoff length.
    // CombinationRule normally uses this.
    // Note that, pair-parameter and per-particle parameter can differ from
    // each other. Lorentz-Bertherot uses the same parameter_type because it is
    // for L-J and L-J-like potentials that has {sigma, epsilon} for each
    // particle and also for each pair of particles.
    template<typename InputIterator>
    real_type max_cutoff(const InputIterator first, const InputIterator last) const noexcept
    {
        static_assert(std::is_same<
                typename std::iterator_traits<InputIterator>::value_type,
                parameter_type>::value, "");

        if(first == last) {return 1;}

        real_type max_sigma = 0;
        for(auto iter = first; iter != last; ++iter)
        {
            const auto& parameter = *iter;
            max_sigma = std::max(max_sigma, parameter.sigma);
        }
        return max_sigma * cutoff_ratio_;
    }
    // It returns absolute cutoff length using pair-parameter.
    // `CombinationTable` uses this.
    real_type absolute_cutoff(const parameter_type& params) const noexcept
    {
        return params.sigma * cutoff_ratio_;
    }

    real_type cutoff_ratio() const noexcept {return cutoff_ratio_;}

  private:

    real_type cutoff_ratio_;
    real_type coef_at_cutoff_;
};

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class LennardJonesAttractivePotential<double>;
extern template class LennardJonesAttractivePotential<float >;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
