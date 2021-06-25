#ifndef MJOLNIR_FORCEFIELD_GLOBAL_UNIFORM_LENNARD_JONES_POTENTIAL_HPP
#define MJOLNIR_FORCEFIELD_GLOBAL_UNIFORM_LENNARD_JONES_POTENTIAL_HPP
#include <mjolnir/util/empty.hpp>
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
class UniformLennardJonesPotential
{
  public:
    using real_type = realT;
    using parameter_type = empty_t;

    static constexpr real_type default_cutoff() noexcept
    {
        return real_type(2.5);
    }

  public:

    explicit UniformLennardJonesPotential(const real_type& cutoff_ratio,
            const real_type sigma, const real_type epsilon) noexcept
        : cutoff_ratio_(cutoff_ratio), coef_at_cutoff_(
                std::pow(real_type(1) / cutoff_ratio, 12) -
                std::pow(real_type(1) / cutoff_ratio,  6)),
          abs_cutoff_(sigma * cutoff_ratio), sigma_(sigma),
          epsilon_(epsilon), epsilon4_(epsilon * 4), epsilon24_(epsilon * 24)
    {}
    ~UniformLennardJonesPotential() = default;

    real_type potential(const real_type r, const parameter_type&) const noexcept
    {
        if(this->abs_cutoff_ < r){return real_type(0);}

        const real_type sr1 = this->sigma_ / r;
        const real_type sr3 = sr1 * sr1 * sr1;
        const real_type sr6 = sr3 * sr3;

        return this->epsilon4_ * (sr6 * (sr6 - real_type(1)) - this->coef_at_cutoff_);
    }
    real_type derivative(const real_type r, const parameter_type&) const noexcept
    {
        if(this->abs_cutoff_ < r){return real_type(0);}

        const real_type rinv = 1 / r;
        const real_type sr1  = this->sigma_ * rinv;
        const real_type sr3  = sr1 * sr1 * sr1;
        const real_type sr6  = sr3 * sr3;

        return this->epsilon24_ * sr6 * (real_type(1) - real_type(2) * sr6) * rinv;
    }

    template<typename T>
    void initialize(const System<T>&) noexcept {return;}

    template<typename T>
    void update(const System<T>&) noexcept {return;}

    real_type cutoff_ratio()   const noexcept {return cutoff_ratio_;}
    real_type coef_at_cutoff() const noexcept {return coef_at_cutoff_;}

    real_type sigma()   const noexcept {return sigma_;}
    real_type epsilon() const noexcept {return epsilon_;}

    static const char* name() noexcept {return "LennardJones";}

    // It takes per-particle parameters and return the maximum cutoff length.
    // CombinationRule normally uses this.
    // Note that, pair-parameter and per-particle parameter can differ from
    // each other. Lorentz-Bertherot uses the same parameter_type because it is
    // for L-J and L-J-like potentials that has {sigma, epsilon} for each
    // particle and also for each pair of particles.
    template<typename InputIterator>
    real_type max_cutoff(const InputIterator, const InputIterator) const noexcept
    {
        return this->abs_cutoff_;
    }
    // It returns absolute cutoff length using pair-parameter.
    // `CombinationTable` uses this.
    real_type absolute_cutoff(const parameter_type&) const noexcept
    {
        return this->abs_cutoff_;
    }

  private:

    real_type cutoff_ratio_;
    real_type coef_at_cutoff_;
    real_type abs_cutoff_;
    real_type sigma_;
    real_type epsilon_;
    real_type epsilon4_;
    real_type epsilon24_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class UniformLennardJonesPotential<double>;
extern template class UniformLennardJonesPotential<float >;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
