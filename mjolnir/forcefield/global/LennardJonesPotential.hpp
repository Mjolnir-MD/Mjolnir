#ifndef MJOLNIR_POTENTIAL_GLOBAL_LENNARD_JONES_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_LENNARD_JONES_POTENTIAL_HPP
#include <algorithm>
#include <utility>
#include <cmath>

namespace mjolnir
{

template<typename T> class System;

// Well-known Lennard-Jones interaction with Lorentz-Berthelot combining rules.
// This class contains sigmas and epsilons of the particles and calculates
// energy and derivative of the potential function.
template<typename realT>
class LennardJonesPotential
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

    explicit LennardJonesPotential(const real_type& cutoff_ratio) noexcept
        : cutoff_ratio_(cutoff_ratio), coef_at_cutoff_(
                std::pow(real_type(1) / cutoff_ratio, 12) -
                std::pow(real_type(1) / cutoff_ratio,  6))
    {}
    ~LennardJonesPotential() = default;

    real_type potential(const real_type r, const parameter_type& params) const noexcept
    {
        if(params.sigma * this->cutoff_ratio_ < r){return real_type(0);}

        const real_type sr1 = params.sigma / r;
        const real_type sr3 = sr1 * sr1 * sr1;
        const real_type sr6 = sr3 * sr3;

        return real_type(4) * params.epsilon *
                (sr6 * (sr6 - real_type(1)) - this->coef_at_cutoff_);
    }
    real_type derivative(const real_type r, const parameter_type& params) const noexcept
    {
        if(params.sigma * this->cutoff_ratio_ < r){return real_type(0);}

        const real_type rinv = 1 / r;
        const real_type sr1  = params.sigma * rinv;
        const real_type sr3  = sr1 * sr1 * sr1;
        const real_type sr6  = sr3 * sr3;

        return real_type(24) * params.epsilon *
                sr6 * (real_type(1) - real_type(2) * sr6) * rinv;
    }

    template<typename T>
    void initialize(const System<T>&) noexcept {return;}

    template<typename T>
    void update(const System<T>&) noexcept {return;}

    real_type cutoff_ratio()   const noexcept {return cutoff_ratio_;}
    real_type coef_at_cutoff() const noexcept {return coef_at_cutoff_;}

    static const char* name() noexcept {return "LennardJones";}

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

  private:

    real_type cutoff_ratio_;
    real_type coef_at_cutoff_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class LennardJonesPotential<double>;
extern template class LennardJonesPotential<float >;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
