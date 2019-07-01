#ifndef MJOLNIR_POTENTIAL_LOCAL_FLEXIBLE_LOCAL_ANGLE_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_LOCAL_FLEXIBLE_LOCAL_ANGLE_POTENTIAL_HPP
#include <array>
#include <algorithm>
#include <limits>
#include <cassert>
#include <cmath>

namespace mjolnir
{

template<typename T> class System;

// Flexible Local Angle potential (T. Terakawa and S. Takada Biophys J 2011)
// NOTE: It assumes that each theta value in histogram is same as the default.
//     : It requires the unit to be [kcal/mol] & [angstrom]!
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
    // TODO: support unit systems
    static constexpr real_type max_force  =  30.0;
    static constexpr real_type min_force  = -30.0;

    static constexpr std::array<real_type, 10> thetas = {
        {1.30900, 1.48353, 1.65806, 1.83260, 2.00713,
         2.18166, 2.35619, 2.53073, 2.70526, 2.87979}
    };
    static constexpr real_type dtheta  = (2.87979 - 1.30900) / 9.0;
    static constexpr real_type rdtheta = 1.0 / dtheta;

  public:

    FlexibleLocalAnglePotential(const real_type k,
                                const std::array<real_type, 10>& ys,
                                const std::array<real_type, 10>& d2ys)
        : min_theta(thetas.front()), max_theta(thetas.back()),
          k_(k), ys_(ys), d2ys_(d2ys)
    {
        // set the range and the parameters from table
        // from cafemol3/mloop_flexible_local.F90
        real_type th = thetas[0];
        const real_type center_th = (max_theta + min_theta) * 0.5;
        this->min_theta_ene = spline_interpolate(min_theta);
        this->max_theta_ene = spline_interpolate(max_theta - 1e-4);

        this->min_energy = min_theta_ene;
        while(th < max_theta)
        {
            const real_type energy = this->spline_interpolate(th);
            const real_type force  = this->spline_derivative(th);
            this->min_energy = std::min(min_energy, energy);

            if(force < min_force)
            {
                this->min_theta     = th;
                this->min_theta_ene = energy;
            }
            if(max_force < force && center_th < th && max_theta == thetas.back())
            {
                this->max_theta     = th;
                this->max_theta_ene = energy;
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
            return k_ * (this->spline_interpolate(th) - min_energy);
        }
    }

    real_type derivative(const real_type th) const noexcept
    {
        if     (th <  min_theta) {return min_force;}
        else if(th >= max_theta) {return max_force;}
        else {return this->spline_derivative(th) * k_;}
    }

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() noexcept {return "FlexibleLocalAngle";}

    real_type                        k()   const noexcept {return k_;}
    std::array<real_type, 10> const& y()   const noexcept {return ys_;}
    std::array<real_type, 10> const& d2y() const noexcept {return d2ys_;}

    real_type cutoff() const noexcept // no cutoff exists.
    {return std::numeric_limits<real_type>::infinity();}

  private:

    real_type spline_interpolate(const real_type th) const noexcept
    {
        constexpr real_type one_over_six = real_type(1.0) / real_type(6.0);

        std::size_t n = std::floor((th - min_theta) * rdtheta);

        if     (thetas[n+1] < th) {n++;}
        else if(th <   thetas[n]) {n--;}
        assert(n < 9);
        assert(thetas[n] <= th && th <= thetas[n+1]);

        const real_type a = (thetas[n+1] - th) * rdtheta;
        const real_type b = (th - thetas[n  ]) * rdtheta;

        const real_type e1 = a * ys_[n] + b * ys_[n+1];
        const real_type e2 =
            ((a * a * a - a) * d2ys_[n] + (b * b * b - b) * d2ys_[n+1]) *
            dtheta * dtheta * one_over_six;

        return e1 + e2;
    }

    real_type spline_derivative(const real_type th) const noexcept
    {
        constexpr real_type one_over_six = real_type(1.0) / real_type(6.0);

        std::size_t n = std::floor((th - min_theta) * rdtheta);

        if     (thetas[n+1] < th) {n++;}
        else if(th <   thetas[n]) {n--;}
        assert(n < 9);
        assert(thetas[n] <= th && th <= thetas[n+1]);

        const real_type a = (thetas[n+1] - th) * rdtheta;
        const real_type b = (th - thetas[n  ]) * rdtheta;

        const real_type f1 = (ys_[n+1] - ys_[n]) * rdtheta;
        const real_type f2 = (
            (3 * b * b - 1) * d2ys_[n+1] - (3 * a * a - 1) * d2ys_[n]
            ) * dtheta * one_over_six;

        return f1 + f2;
    }

  private:
    real_type min_energy;
    real_type min_theta;
    real_type max_theta;
    real_type min_theta_ene;
    real_type max_theta_ene;

    real_type k_;
    std::array<real_type, 10> ys_;
    std::array<real_type, 10> d2ys_;
};
template<typename realT>
constexpr realT FlexibleLocalAnglePotential<realT>::max_force;
template<typename realT>
constexpr realT FlexibleLocalAnglePotential<realT>::min_force;
template<typename realT>
constexpr realT FlexibleLocalAnglePotential<realT>::dtheta;
template<typename realT>
constexpr realT FlexibleLocalAnglePotential<realT>::rdtheta;
template<typename realT>
constexpr std::array<realT, 10> FlexibleLocalAnglePotential<realT>::thetas;

} // mjolnir
#endif // MJOLNIR_FLEXIBLE_LOCAL_ANGLE_POTENTIAL
