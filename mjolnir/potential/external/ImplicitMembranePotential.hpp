#ifndef MJOLNIR_IMPLICIT_MEMBRANE_POTENTIAL
#define MJOLNIR_IMPLICIT_MEMBRANE_POTENTIAL
#include <mjolnir/core/System.hpp>
#include <cmath>

namespace mjolnir
{
/* Implicit membrane potential                                                *
 *  V(z) = k * tanh(bend * (|z| - thick/2))                                   *
 * dV/dr = (z/|z|) * k * (cosh^2(bend * (|z| - thick/2)))                     *
 *                 y                                                          *
 *  _________      ^      _______                                             *
 *  _________\_____|_____/______x                                             *
 *  -thick/2  \____|____/ thick/2                                             *
 *               -k|                                                          *
 * Cutoff ratio ensure 1/1000 accuracy.                                       */
template<typename realT>
class ImplicitMembranePotential
{
  public:
    using real_type = realT;

    static constexpr real_type cutoff_ratio = 4.0;

  public:
    ImplicitMembranePotential(
        const real_type thick, const real_type magnitude, const real_type bend,
        std::vector<real_type> hydrophobicities)
        : half_thick_(thick / 2), interaction_magnitude_(magnitude),
          bend_(bend), hydrophobicities_(std::move(hydrophobicities))
    {}
    ~ImplicitMembranePotential() = default;

    real_type potential(const std::size_t i, const real_type z) const noexcept
    {
        return hydrophobicities_[i] * interaction_magnitude_ *
            std::tanh(bend_ * (std::abs(z) - half_thick_));
    }
    real_type derivative(const std::size_t i, const real_type z) const noexcept
    {
        return hydrophobicities_[i] * std::copysign(real_type(1.0), z) *
            interaction_magnitude_ * bend_ /
            std::pow((std::cosh(bend_ * (std::abs(z) - half_thick_))), 2);
    }
    real_type max_cutoff_length() const noexcept
    {
        return cutoff_ratio / this->bend_ + this->half_thick_;
    }

    // nothing to be done if system parameter (e.g. temperature) changes
    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    std::vector<std::size_t> participants() const
    {
        std::vector<std::size_t> retval;
        retval.reserve(this->hydrophobicities_.size());
        for(std::size_t i=0; i<this->hydrophobicities_.size(); ++i)
        {
            if(this->hydrophobicities_[i] != real_type(0))
            {
                retval.push_back(i);
            }
        }
        return retval;
    }

    const char* name() const noexcept {return "ImplicitMembrane";}

    real_type  half_thick() const noexcept {return half_thick_;}
    real_type& half_thick()       noexcept {return half_thick_;}

    real_type  bend() const noexcept {return bend_;}
    real_type& bend()       noexcept {return bend_;}

    real_type  interaction_magnitude() const noexcept
    {return interaction_magnitude_;}
    real_type& interaction_magnitude()       noexcept
    {return interaction_magnitude_;}

    std::vector<real_type> const&
    hydrophobicities() const noexcept {return hydrophobicities_;}

  private:

    real_type half_thick_;            // half of thickness of the membrane.
    real_type interaction_magnitude_; // overall scaling parameter.
    real_type bend_;                  // the slope of tanh curve.
    std::vector<real_type> hydrophobicities_;
};

template<typename realT>
constexpr typename ImplicitMembranePotential<realT>::real_type
ImplicitMembranePotential<realT>::cutoff_ratio;

} // mjolnir
#endif /* MJOLNIR_IMPLICIT_MEMBRANE_POTENTIAL */
