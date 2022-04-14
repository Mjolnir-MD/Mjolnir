#ifndef MJOLNIR_POTENTIAL_GLOBAL_CUBIC_FUNCTION_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_CUBIC_FUNCTION_POTENTIAL_HPP
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/util/empty.hpp>
#include <mjolnir/util/logger.hpp>
#include <vector>

namespace mjolnir
{

template<typename T> class System;

// This potential is same to GlobalStoichiometricInteractionPotential
// and exist to use that in non-stoichiometric interaction.
template<typename realT>
class UniformCubicFunctionPotential
{
  public:
    using real_type           = realT;
    using parameter_type      = empty_t; // no particle_specific parameter

  public:

    explicit UniformCubicFunctionPotential(const real_type epsilon,
            const real_type v0, const real_type range) noexcept
        : epsilon_(epsilon), v0_(v0), range_(range), v0_range_(v0 + range),
          inv_range_(1.0 / range), inv_range_6_(6.0 / range)
    {}
    ~UniformCubicFunctionPotential() = default;

    // forwarding functions for clarity...
    real_type potential(const real_type r, const parameter_type&) const noexcept
    {
        if(r < v0_){return -1.0;}
        else if(v0_range_ < r){return 0.0;}

        const real_type r_v0           = r - v0_;
        const real_type r_v0_inv_range = r_v0 * inv_range_;
        const real_type r_v0_inv_range_2  = r_v0_inv_range * r_v0_inv_range;
        const real_type r_v0_inv_range_3  = r_v0_inv_range_2 * r_v0_inv_range;

        return -epsilon_ * (1.0 - 3.0 * r_v0_inv_range_2 + 2.0 * r_v0_inv_range_3);
    }
    real_type derivative(const real_type r, const parameter_type&) const noexcept
    {
        if(r < v0_ || v0_range_ < r){return 0.0;}

        const real_type r_v0          = r - v0_;
        const real_type r_v0_inv_range   = r_v0 * inv_range_;
        const real_type r_v0_inv_range_2 = r_v0_inv_range * r_v0_inv_range;

        return -epsilon_ * inv_range_6_ * (r_v0_inv_range_2 - r_v0_inv_range);
    }


    template<typename T>
    void initialize(const System<T>&) noexcept {return;}

    template<typename T>
    void update(const System<T>&) noexcept {return;}

    real_type epsilon()           const noexcept {return epsilon_;}
    real_type v0()                const noexcept {return v0_;}
    real_type interaction_range() const noexcept {return range_;}

    static const char* name() noexcept {return "UniformCubicFunction";}

    // It takes per-particle parameters and return the maximum cutoff length.
    template<typename InputIterator>
    real_type max_cutoff(const InputIterator, const InputIterator) const noexcept
    {
        return this->v0_range_;
    }
    real_type absolute_cutoff(const parameter_type&) const noexcept
    {
        return this->v0_range_;
    }

  private:

    real_type epsilon_;
    real_type v0_;
    real_type range_;
    real_type v0_range_;
    real_type inv_range_;
    real_type inv_range_6_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class UniformCubicFunctionPotential<double>;
extern template class UniformCubicFunctionPotential<float >;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif /* MJOLNIR_POTENTIAL_GLOBAL_CUBIC_FUNCTION_POTENTIAL_HPP */
