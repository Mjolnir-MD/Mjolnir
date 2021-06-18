#ifndef MJOLNIR_POTENTIAL_GLOBAL_WCA_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_WCA_POTENTIAL_HPP
#include <algorithm>
#include <utility>
#include <cmath>

namespace mjolnir
{

template<typename T> class System;

// Well-known WCA interaction with Lorentz-Berthelot combining rules.
// This class contains sigmas and epsilons of the particles and calculates
// energy and derivative of the potential function.
template<typename realT>
class WCAPotential
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
        return real_type(1.12246204831); // pow(2.0, 1.0 / 6.0)
    }

  public:

    WCAPotential() noexcept {}
    ~WCAPotential() = default;

    real_type potential(const real_type r, const parameter_type& params) const noexcept
    {
        if(params.sigma * default_cutoff() < r){return real_type(0);}

        const real_type r1s1 = params.sigma / r;
        const real_type r3s3 = r1s1 * r1s1 * r1s1;
        const real_type r6s6 = r3s3 * r3s3;
        return real_type(4) * params.epsilon * (r6s6 * (r6s6 - real_type(1)) + real_type(0.25));
    }
    real_type derivative(const real_type r, const parameter_type& params) const noexcept
    {
        if(params.sigma * default_cutoff() < r){return real_type(0);}

        const real_type rinv = 1 / r;
        const real_type r1s1 = params.sigma * rinv;
        const real_type r3s3 = r1s1 * r1s1 * r1s1;
        const real_type r6s6 = r3s3 * r3s3;

        return real_type(24) * params.epsilon * (r6s6 * (real_type(1) - 2 * r6s6)) * rinv;
    }

    template<typename T>
    void initialize(const System<T>&) noexcept {return;}

    template<typename T>
    void update(const System<T>&) noexcept {return;}

    real_type cutoff_ratio()   const noexcept {return default_cutoff();}
    real_type coef_at_cutoff() const noexcept {return real_type(0);}

    static const char* name() noexcept {return "WCA";}

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
        return max_sigma * default_cutoff();
    }
    // It returns absolute cutoff length using pair-parameter.
    // `CombinationTable` uses this.
    real_type absolute_cutoff(const parameter_type& params) const noexcept
    {
        return params.sigma * default_cutoff();
    }
};
} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class WCAPotential<double>;
extern template class WCAPotential<float >;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
