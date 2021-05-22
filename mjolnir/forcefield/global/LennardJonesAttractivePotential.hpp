#ifndef MJOLNIR_POTENTIAL_GLOBAL_LENNARD_JONES_ATTRACTIVE_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_LENNARD_JONES_ATTRACTIVE_POTENTIAL_HPP
#include <utility>
#include <cmath>

namespace mjolnir
{

template<typename T> struct System;

template<typename realT>
class LennardJonesAttractivePotential
{
  public:
    using real_type      = realT;
    using parameter_type = std::pair<real_type, real_type>; // {sigma, epsilon}
    using self_type      = LennardJonesAttractivePotential<real_type>;

    static constexpr real_type default_cutoff() noexcept
    {
        return real_type(2.5);
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{real_type(0), real_type(0)};
    }

    static void set_cutoff_ratio(const real_type ratio)
    {
        if(self_type::cutoff_ratio < ratio)
        {
            self_type::cutoff_ratio  = ratio;
            self_type::coef_at_cutoff = std::pow(real_type(1) / ratio, 12)
                                      - std::pow(real_type(1) / ratio, 6);
        }
        return;
    }

    static real_type cutoff_ratio;
    static real_type coef_at_cutoff;

  public:

    LennardJonesAttractivePotential(const real_type sigma, const real_type epsilon) noexcept
        : sigma_(sigma), epsilon_(epsilon)
    {}
    explicit LennardJonesAttractivePotential(const parameter_type& params) noexcept
        : sigma_(params.first), epsilon_(params.second)
    {}
    ~LennardJonesAttractivePotential() = default;

    real_type potential(const real_type r) const noexcept
    {
        constexpr real_type sixth_root_of_two(1.12246204831);

        if(r < this->sigma_ * sixth_root_of_two)
        {
            return -epsilon_;
        }
        else if(this->sigma_ * self_type::cutoff_ratio < r)
        {
            return 0;
        }

        const real_type sr1 = this->sigma_ / r;
        const real_type sr3 = sr1 * sr1 * sr1;
        const real_type sr6 = sr3 * sr3;
        return 4 * this->epsilon_ * (sr6 * (sr6 - 1) - self_type::coef_at_cutoff);
    }
    real_type derivative(const real_type r) const noexcept
    {
        constexpr real_type sixth_root_of_two(1.12246204831);

        if(r < sigma_ * sixth_root_of_two ||
               sigma_ * self_type::cutoff_ratio < r)
        {
            return 0;
        }

        const real_type rinv = 1 / r;
        const real_type sr1 = sigma_ * rinv;
        const real_type sr3 = sr1 * sr1 * sr1;
        const real_type sr6 = sr3 * sr3;
        return 24 * epsilon_ * (sr6 - 2 * sr6 * sr6) * rinv;
    }


    template<typename T>
    void initialize(const System<T>&) noexcept {return;}

    template<typename T>
    void update(const System<T>&) noexcept {return;}

    static const char* name() noexcept {return "LennardJonesAttractive";}

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
typename LennardJonesAttractivePotential<realT>::real_type LennardJonesAttractivePotential<realT>::cutoff_ratio  = 0.0;
template<typename realT>
typename LennardJonesAttractivePotential<realT>::real_type LennardJonesAttractivePotential<realT>::coef_at_cutoff = 0.0;

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
