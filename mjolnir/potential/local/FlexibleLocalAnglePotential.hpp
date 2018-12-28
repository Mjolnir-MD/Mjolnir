#ifndef MJOLNIR_FLEXIBLE_LOCAL_ANGLE_POTENTIAL
#define MJOLNIR_FLEXIBLE_LOCAL_ANGLE_POTENTIAL
#include <array>
#include <algorithm>
#include <cassert>
#include <cmath>

namespace mjolnir
{

template<typename T> class System;

// Flexible Local Angle potential (T. Terakawa and S. Takada Biophys J 2011)
// NOTE: It assumes that each theta value in histogram is same as the default.
//
// Here, the implementation is derived from the original implementation in
// the CGMD package, Cafemol (Kenzaki et al., JCTC 2011).
// As the original implementation does, the first term represents the value of
// energy, and the second term represents the value of the *second* derivative
// of the energy.
template<typename realT>
class FlexibleLocalAnglePotential
{
  public:
    using real_type = realT;
    static constexpr real_type max_force  =  30.0;
    static constexpr real_type min_force  = -30.0;

  public:
    FlexibleLocalAnglePotential(const real_type k,
            const std::array<real_type, 10>& term1,
            const std::array<real_type, 10>& term2)
        : min_theta(1.30900), max_theta(2.87979),
          dtheta((max_theta - min_theta) / 9.0), inv_dtheta(1. / dtheta), k_(k),
          thetas{{1.30900, 1.48353, 1.65806, 1.83260, 2.00713,
                  2.18166, 2.35619, 2.53073, 2.70526, 2.87979}},
          term1_(term1), term2_(term2)
    {
        // from cafemol3/mloop_flexible_local.F90
        real_type th = thetas[0];
        const real_type center_th  = (max_theta + min_theta) * 0.5;
        min_theta_ene = min_energy = spline_interpolate(min_theta);
        max_theta_ene =              spline_interpolate(max_theta - 1e-4);

        while(th < max_theta)
        {
            const real_type energy = spline_interpolate(th);
            const real_type force  = spline_derivative(th);

            min_energy = std::min(min_energy, energy);

            if(force < min_force)
            {
                min_theta     = th;
                min_theta_ene = energy;
            }
            if(max_force < force && center_th < th && max_theta == thetas[9])
            {
                max_theta     = th;
                max_theta_ene = energy;
            }
            th += 1e-4;
        }
    }

    ~FlexibleLocalAnglePotential() = default;

    real_type potential(const real_type th) const noexcept
    {
        if(th < min_theta)
        {
            return ((min_force * th + min_theta_ene - min_force * min_theta) -
                     min_energy) * k_;
        }
        else if(th >= max_theta)
        {
            return ((max_force * th + max_theta_ene - max_force * max_theta) -
                     min_energy) * k_;
        }
        else
        {
            return k_ * (spline_interpolate(th) - min_energy);
        }
    }

    real_type derivative(const real_type th) const noexcept
    {
             if(th <  min_theta) {return min_force;}
        else if(th >= max_theta) {return max_force;}
        else {return spline_derivative(th) * k_;}
    }

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() noexcept {return "FlexibleLocalAngle";}

    real_type                        k()   const noexcept {return k_;}
    std::array<real_type, 10> const& y()   const noexcept {return term1_;}
    std::array<real_type, 10> const& d2y() const noexcept {return term2_;}

  private:

    real_type spline_interpolate(const real_type th) const noexcept
    {
        const std::size_t n = std::floor((th - min_theta) * inv_dtheta);
        assert(n < 9);
        const real_type a = (thetas[n+1] - th) * inv_dtheta;
        const real_type b = (th - thetas[n  ]) * inv_dtheta;

        const real_type e1 = a * term1_[n] + b * term1_[n+1];
        const real_type e2 =
            ((a * a * a - a) * term2_[n] + (b * b * b - b) * term2_[n+1]) *
            dtheta * dtheta / real_type(6.0);

        return e1 + e2;
    }

    real_type spline_derivative(const real_type th) const noexcept
    {
        const std::size_t n = std::floor((th - min_theta) * inv_dtheta);
        assert(n < 9);
        const real_type a = (thetas[n+1] - th) * inv_dtheta;
        const real_type b = (th - thetas[n  ]) * inv_dtheta;

        const real_type f1 = (term1_[n+1] - term1_[n]) * inv_dtheta;
        const real_type f2 = (
            (3 * b * b - 1) * term2_[n+1] - (3 * a * a - 1) * term2_[n]
            ) * dtheta / real_type(6.0);

        return f1 + f2;
    }

  private:
    real_type min_energy;
    real_type min_theta_ene;
    real_type max_theta_ene;
    real_type min_theta;
    real_type max_theta;

    real_type dtheta;
    real_type inv_dtheta;
    real_type k_;
    std::array<real_type, 10> thetas;
    std::array<real_type, 10> term1_;
    std::array<real_type, 10> term2_;
};

} // mjolnir
#endif // MJOLNIR_FLEXIBLE_LOCAL_ANGLE_POTENTIAL
