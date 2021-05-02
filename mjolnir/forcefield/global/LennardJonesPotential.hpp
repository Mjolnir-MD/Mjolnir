#ifndef MJOLNIR_POTENTIAL_GLOBAL_LENNARD_JONES_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_LENNARD_JONES_POTENTIAL_HPP
#include <limits>
#include <utility>
#include <cmath>

namespace mjolnir
{

template<typename T> struct System;

// Well-known Lennard-Jones interaction with Lorentz-Berthelot combining rules.
// This class contains sigmas and epsilons of the particles and calculates
// energy and derivative of the potential function.
template<typename realT>
class LennardJonesPotential
{
  public:
    using real_type      = realT;
    using parameter_type = std::pair<real_type, real_type>; // {sigma, epsilon}
    using self_type      = LennardJonesPotential<realT>;

    static constexpr real_type default_cutoff() noexcept
    {
        return real_type(2.5);
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{real_type(0), real_type(0)};
    }

    static void set_cutoff_length(const real_type cutoff_ratio)
    {
        if(self_type::cutoff_length < cutoff_ratio)
        {
            self_type::cutoff_length  = cutoff_ratio;
            self_type::coef_at_cutoff = std::pow(real_type(1) / cutoff_ratio, 12)
                                      - std::pow(real_type(1) / cutoff_ratio, 6);
        }
        return;
    }

    static real_type cutoff_length;
    static real_type coef_at_cutoff;

  public:

    LennardJonesPotential(const real_type sigma, const real_type epsilon) noexcept
        : sigma_(sigma), epsilon_(epsilon)
    {}
    explicit LennardJonesPotential(const parameter_type& params) noexcept
        : sigma_(params.first), epsilon_(params.second)
    {}
    ~LennardJonesPotential() = default;

    real_type potential(const real_type r) const noexcept
    {
        if(this->sigma_ * self_type::cutoff_length < r){return real_type(0);}

        const real_type sr1 = sigma_ / r;
        const real_type sr3 = sr1 * sr1 * sr1;
        const real_type sr6 = sr3 * sr3;

        return real_type(4) * epsilon_ * (sr6 * (sr6 - real_type(1)) - self_type::coef_at_cutoff);
    }
    real_type derivative(const real_type r) const noexcept
    {
        if(this->sigma_ * self_type::cutoff_length < r){return real_type(0);}

        const real_type rinv = 1 / r;
        const real_type sr1  = sigma_ * rinv;
        const real_type sr3  = sr1 * sr1 * sr1;
        const real_type sr6  = sr3 * sr3;

        return real_type(24) * epsilon_ * sr6 * (real_type(1) - real_type(2) * sr6) * rinv;
    }

    template<typename T>
    void initialize(const System<T>&) noexcept {return;}

    template<typename T>
    void update(const System<T>&) noexcept {return;}

    static const char* name() noexcept {return "LennardJones";}

    real_type sigma()   const noexcept {return this->sigma_;}
    real_type epsilon() const noexcept {return this->epsilon_;}

    real_type cutoff()  const noexcept
    {
        return this->sigma_ * self_type::cutoff_length;
    }

  public:

    // To culculate cutoff distance, we need to find the maximum sigma in the
    // existing parameters. But the list of parameters will be given in a variety
    // of ways, like Lorentz-Bertherot rule, combination table, or another way
    // of combination rules.
    //     To find the maximum parameter, we need to provide a way to compare
    // parameters. But the way depends on the functional form of a potential.
    // So this comparator should be defined in a Potential class.
    struct parameter_comparator
    {
        constexpr bool
        operator()(const parameter_type& lhs, const parameter_type& rhs) const noexcept
        {
            return lhs.first < rhs.first; // take larger sigma
        }
    };

    parameter_comparator get_parameter_comparator() const noexcept
    {
        return parameter_comparator{};
    }

  private:

    real_type sigma_;
    real_type epsilon_;
};

template<typename realT>
typename LennardJonesPotential<realT>::real_type LennardJonesPotential<realT>::cutoff_length  = 0.0;
template<typename realT>
typename LennardJonesPotential<realT>::real_type LennardJonesPotential<realT>::coef_at_cutoff = 0.0;

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class LennardJonesPotential<double>;
extern template class LennardJonesPotential<float >;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
