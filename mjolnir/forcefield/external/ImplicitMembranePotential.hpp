#ifndef MJOLNIR_POTENTIAL_EXTERNAL_IMPLICIT_MEMBRANE_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_EXTERNAL_IMPLICIT_MEMBRANE_POTENTIAL_HPP
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
    using parameter_type = real_type;
    using container_type = std::vector<parameter_type>;

    static constexpr real_type default_cutoff() noexcept
    {
        return real_type(4);
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return real_type(0);
    }

  public:

    ImplicitMembranePotential(const real_type thick, const real_type k,
        const real_type bend, const real_type cutoff_ratio,
        const std::vector<std::pair<std::size_t, real_type>>& params)
        : half_thick_(thick / 2), k_(k), bend_(bend), cutoff_ratio_(cutoff_ratio)
    {
        this->parameters_.resize(params.size());
        this->participants_.reserve(params.size());
        for(const auto& idxp: params)
        {
            const auto idx = idxp.first;
            this->participants_.push_back(idx);
            if(idx >= this->parameters_.size())
            {
                this->parameters_.resize(idx+1, default_parameter());
            }
            this->parameters_.at(idx) = idxp.second;
        }
    }
    ~ImplicitMembranePotential() = default;

    real_type potential(const std::size_t i, const real_type z) const noexcept
    {
        return parameters_[i] * k_ * std::tanh(bend_ * (std::abs(z) - half_thick_));
    }
    real_type derivative(const std::size_t i, const real_type z) const noexcept
    {
        return parameters_[i] * std::copysign(real_type(1.0), z) * k_ * bend_ /
            std::pow(std::cosh(bend_ * (std::abs(z) - half_thick_)), 2);
    }
    real_type max_cutoff_length() const noexcept
    {
        return this->cutoff_ratio_ / this->bend_ + this->half_thick_;
    }

    // nothing to be done if system parameter (e.g. temperature) changes
    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    std::vector<std::size_t> const& participants() const noexcept
    {
        return participants_;
    }

    const char* name() const noexcept {return "ImplicitMembrane";}

    real_type  half_thick() const noexcept {return half_thick_;}
    real_type& half_thick()       noexcept {return half_thick_;}

    real_type  bend() const noexcept {return bend_;}
    real_type& bend()       noexcept {return bend_;}

    real_type  k() const noexcept {return k_;}
    real_type& k()       noexcept {return k_;}

    real_type  cutoff() const noexcept {return cutoff_ratio_;}
    real_type& cutoff()       noexcept {return cutoff_ratio_;}

    container_type const& parameters() const noexcept {return parameters_;}
    container_type&       parameters()       noexcept {return parameters_;}

  private:

    real_type half_thick_; // half of thickness of the membrane.
    real_type k_;          // overall scaling parameter.
    real_type bend_;       // the slope of tanh curve.
    real_type cutoff_ratio_;
    container_type           parameters_;
    std::vector<std::size_t> participants_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class ImplicitMembranePotential<double>;
extern template class ImplicitMembranePotential<float>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif /* MJOLNIR_IMPLICIT_MEMBRANE_POTENTIAL */
