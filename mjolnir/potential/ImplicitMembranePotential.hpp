#ifndef MJOLNIR_IMPLICIT_MEMBRANE_POTENTIAL
#define MJOLNIR_IMPLICIT_MEMBRANE_POTENTIAL
#include <cmath>

namespace mjolnir
{
  
  /* Implicit membrane potential & derivative      *
   * potential field dependent on z coordinate.    *
   * tanh is used to represent membrane potential. *
   *  V(z) = ma * tanh(be * (|z| - th))            *
   * dV/dr = ma * z / (|z| * tanh(be * (|z| - th)))*///TODO This is incorrect. fix!
  template<typename traitT>
  class ImplicitMembranePotential
  {
  public:
    typedef traitT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

    //TODO set cutoff_ratio;
    constexpr static real_type cutoff_ratio
    
  public:
    ImplicitMembranePotential(const real_type th, const real_type ma, const real_type be)
      : thickness_(th),interaction_magnitude_(ma),bend_(be);
    {}
    ~ImplicitMembranePotential() = default;

    real_type potential(const real_type z) const
    {
      return this->interaction_magnitude_ * std::tanh(bend_ * (std::abs(z) - thickness_));
    }

    real_type derivative(const real_type z) const
    {
      return this->interaction_magnitude_ * z
	/ (std::abs(z) * std::tanh(bend_ * (std::abs(z) - thickness_)));
    }

    void reset_parameter(const std::string&, const real_type){return;}
    
  private:
    
    const real_type thickness_;//membrane thickness.
    const real_type interaction_magnitude_;
    const real_type bend_;//bend_ decide the slope of tanh carve.
  };
}
#endif /* MJOLNIR_IMPLICIT_MEMBRANE_POTENTIAL */
