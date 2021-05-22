#ifndef MJOLNIR_POTENTIAL_GLOBAL_WCA_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_WCA_POTENTIAL_HPP
#include <utility>
#include <cmath>

namespace mjolnir
{

template<typename T> struct System;

// Well-known WCA interaction with Lorentz-Berthelot combining rules.
// This class contains sigmas and epsilons of the particles and calculates
// energy and derivative of the potential function.
template<typename realT>
class WCAPotential
{
  public:
    using real_type      = realT;
    using parameter_type = std::pair<real_type, real_type>; // {sigma, epsilon}
    using self_type      = WCAPotential<real_type>;

    static constexpr real_type default_cutoff() noexcept
    {
        return real_type(1.12246204831); // pow(2.0, 1.0 / 6.0)
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{real_type(0), real_type(0)};
    }

    static void set_cutoff_ratio(const real_type) // does not have any effect.
    {
        return;
    }

    static constexpr real_type cutoff_ratio = real_type(1.12246204831);
    static constexpr real_type coef_at_cutoff = real_type(0);

  public:

    explicit WCAPotential(const parameter_type& params) noexcept
        : sigma_(params.first), epsilon_(params.second)
    {}
    ~WCAPotential() = default;

    real_type potential(const real_type r) const noexcept
    {
        if(this->sigma_ * self_type::cutoff_ratio < r){return real_type(0);}

        const real_type r1s1 = sigma_ / r;
        const real_type r3s3 = r1s1 * r1s1 * r1s1;
        const real_type r6s6 = r3s3 * r3s3;
        return real_type(4) * this->epsilon_ * (r6s6 * (r6s6 - real_type(1)) + real_type(0.25));
    }
    real_type derivative(const real_type r) const noexcept
    {
        if(this->sigma_ * self_type::cutoff_ratio < r){return real_type(0);}

        const real_type rinv = 1 / r;
        const real_type r1s1 = sigma_ * rinv;
        const real_type r3s3 = r1s1 * r1s1 * r1s1;
        const real_type r6s6 = r3s3 * r3s3;

        return real_type(24) * this->epsilon_ * (r6s6 * (real_type(1) - 2 * r6s6)) * rinv;
    }

    template<typename T>
    void initialize(const System<T>&) noexcept {return;}

    template<typename T>
    void update(const System<T>&) noexcept {return;}

    static const char* name() noexcept {return "WCA";}

    real_type sigma()   const noexcept {return this->sigma_;}
    real_type epsilon() const noexcept {return this->epsilon_;}

    real_type cutoff()  const noexcept
    {
        return this->sigma_ * self_type::cutoff_ratio;
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

  private:

    real_type sigma_;
    real_type epsilon_;
};
template<typename realT>
constexpr typename WCAPotential<realT>::real_type WCAPotential<realT>::cutoff_ratio;
template<typename realT>
constexpr typename WCAPotential<realT>::real_type WCAPotential<realT>::coef_at_cutoff;


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
