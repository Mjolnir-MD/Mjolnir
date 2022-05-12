#ifndef MJOLNIR_FORCEFIELD_GLOBAL_STOICHIOMETRIC_INTERACTION_POTENTIAL_HPP
#define MJOLNIR_FORCEFIELD_GLOBAL_STOICHIOMETRIC_INTERACTION_POTENTIAL_HPP
#include <mjolnir/util/empty.hpp>
#include <algorithm>

namespace mjolnir
{

template<typename T> class System;

template<typename realT>
class GlobalStoichiometricInteractionPotential
{
  public:

    using real_type      = realT;
    using parameter_type = empty_t;

  public:

    explicit GlobalStoichiometricInteractionPotential(
            const real_type v0, const real_type range) noexcept
        : v0_(v0), v0_range_(v0 + range),
          inv_range_(1.0 / range), inv_range_6_(6.0 / range)
    {}
    ~GlobalStoichiometricInteractionPotential() = default;
    GlobalStoichiometricInteractionPotential(const GlobalStoichiometricInteractionPotential&) = default;
    GlobalStoichiometricInteractionPotential(GlobalStoichiometricInteractionPotential&&) = default;
    GlobalStoichiometricInteractionPotential& operator=(const GlobalStoichiometricInteractionPotential&) = default;
    GlobalStoichiometricInteractionPotential& operator=(GlobalStoichiometricInteractionPotential&&) = default;

    real_type potential(const real_type r) const noexcept
    {
        if(r < v0_){return 1.0;}
        else if(v0_range_ < r){return 0.0;}

        const real_type r_v0              = r - v0_;
        const real_type r_v0_inv_range    = r_v0 * inv_range_;
        const real_type r_v0_inv_range_2  = r_v0_inv_range * r_v0_inv_range;
        const real_type r_v0_inv_range_3  = r_v0_inv_range_2 * r_v0_inv_range;
        return 1.0 - 3.0 * r_v0_inv_range_2 + 2.0 * r_v0_inv_range_3;
    }

    real_type derivative(const real_type r) const noexcept
    {
        if(r < v0_ || v0_range_ < r){return 0.0;}

        const real_type r_v0             = r - v0_;
        const real_type r_v0_inv_range   = r_v0 * inv_range_;
        const real_type r_v0_inv_range_2 = r_v0_inv_range * r_v0_inv_range;
        return inv_range_6_ * (r_v0_inv_range_2 - r_v0_inv_range);
    }

    real_type max_cutoff_length() const noexcept
    {
        return v0_range_;
    }

    template<typename T>
    void initialize(const System<T>&) noexcept {return;}

    template<typename T>
    void update(const System<T>) noexcept {return;}

    // -----------------------------------------------------------------------
    // used by Observer.
    static const char* name() noexcept {return "GlobalStoichiometricInteraction";}

  private:

    real_type                v0_;
    real_type                v0_range_;
    real_type                inv_range_;
    real_type                inv_range_6_;
};

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize major use-cases
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

namespace mjolnir
{
extern template class GlobalStoichiometricInteractionPotential<double>;
extern template class GlobalStoichiometricInteractionPotential<float >;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif // MJOLNIR_POTENTIAL_GLOBAL_STOICHIOMETRIC_INTERACTION_POTENTIAL_HPP
