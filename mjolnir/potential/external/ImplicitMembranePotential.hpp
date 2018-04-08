#ifndef MJOLNIR_IMPLICIT_MEMBRANE_POTENTIAL
#define MJOLNIR_IMPLICIT_MEMBRANE_POTENTIAL
#include <mjolnir/core/System.hpp>
#include <cmath>

namespace mjolnir
{


/* Implicit membrane potential & derivative           *
 * potential field dependent on z coordinate.         *
 * tanh is used to represent membrane potential.      *
 *  V(z) = ma * tanh(be * (|z| - th/2))               *
 * dV/dr = (z/|z|) * ma * (cosh^2(be * (|z| - th/2))) *
 * Cutoff ratio ensure 1/1000 accuracy.               */
template<typename traitsT>
class ImplicitMembranePotential
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef real_type parameter_type;

    constexpr static real_type cutoff_ratio = 4.0;

  public:
    ImplicitMembranePotential() = default;
    ImplicitMembranePotential(const real_type th, const real_type ma,
        const real_type be, const std::vector<parameter_type>& hydrophobicities)
        : half_thick_(th * 0.5), interaction_magnitude_(ma),
          bend_(be), hydrophobicities_(hydrophobicities)
    {}

    ImplicitMembranePotential(real_type th, real_type ma,
        real_type be, std::vector<parameter_type>&& hydrophobicities)
        : half_thick_(th * 0.5), interaction_magnitude_(ma), bend_(be),
          hydrophobicities_(std::move(hydrophobicities))
    {}
    ~ImplicitMembranePotential() = default;

    real_type  half_thick() const noexcept {return half_thick_;}
    real_type& half_thick()       noexcept {return half_thick_;}

    real_type  interaction_magnitude() const noexcept
    {return interaction_magnitude_;}
    real_type& interaction_magnitude() noexcept
    {return interaction_magnitude_;}

    real_type  bend() const noexcept {return bend_;}
    real_type& bend()       noexcept {return bend_;}

    void hydrophobicities_emplace(const parameter_type);

    void set_hydrophobicities(const std::vector<parameter_type>&);
    void set_hydrophobicities(std::vector<parameter_type>&&);

    std::size_t size() const noexcept {return hydrophobicities_.size();}
    void resize (const std::size_t i){hydrophobicities_.resize(i);}
    void reserve(const std::size_t i){hydrophobicities_.reserve(i);}
    void clear() {hydrophobicities_.clear();}

    parameter_type&       operator[](const std::size_t i) noexcept
    {return hydrophobicities_[i];}
    parameter_type const& operator[](const std::size_t i) const noexcept
    {return hydrophobicities_[i];}
    parameter_type&       at(const std::size_t i)
    {return hydrophobicities_.at(i);}
    parameter_type const& at(const std::size_t i) const
    {return hydrophobicities_.at(i);}

    real_type potential (const std::size_t i, const real_type z) const;
    real_type derivative(const std::size_t i, const real_type z) const;
    real_type max_cutoff_length() const noexcept;

    // nothing to be done if system parameter (e.g. temperature) changes
    void update(const System<traitsT>&) const noexcept {return;}

  private:

    real_type half_thick_;//membrane half of thickness.
    real_type interaction_magnitude_;
    real_type bend_;//bend_ decide the slope of tanh carve.
    std::vector<parameter_type> hydrophobicities_;
};

template<typename traitsT>
inline void
ImplicitMembranePotential<traitsT>::hydrophobicities_emplace(
    const parameter_type hydrophobicity)
{
    hydrophobicities_.emplace_back(hydrophobicity);
    return;
}

template<typename traitsT>
inline void
ImplicitMembranePotential<traitsT>::set_hydrophobicities(
    const std::vector<parameter_type>& hydrophobicities)
{
    hydrophobicities_ = hydrophobicities;
    return;
}

template<typename traitsT>
inline void
ImplicitMembranePotential<traitsT>::set_hydrophobicities(
    std::vector<parameter_type>&& hydrophobicities)
{
    hydrophobicities_ = std::move(hydrophobicities);
    return;
}

template<typename traitsT>
inline typename ImplicitMembranePotential<traitsT>::real_type
ImplicitMembranePotential<traitsT>::potential(
    const std::size_t i, const real_type z) const
{
    return hydrophobicities_[i] * interaction_magnitude_ *
        std::tanh(bend_ * (std::abs(z) - half_thick_));
}

template<typename traitsT>
inline typename ImplicitMembranePotential<traitsT>::real_type
ImplicitMembranePotential<traitsT>::derivative(
    const std::size_t i, const real_type z) const
{
    return hydrophobicities_[i] * std::copysign(1.0, z) *
        interaction_magnitude_ * bend_ /
        std::pow((std::cosh(bend_ * (std::abs(z) - half_thick_))), 2);
}

template<typename traitsT>
inline typename ImplicitMembranePotential<traitsT>::real_type
ImplicitMembranePotential<traitsT>::max_cutoff_length() const noexcept
{
    return cutoff_ratio / bend_ + half_thick_;
}

}
#endif /* MJOLNIR_IMPLICIT_MEMBRANE_POTENTIAL */
