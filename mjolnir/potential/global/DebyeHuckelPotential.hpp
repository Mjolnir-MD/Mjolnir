#ifndef MJOLNIR_DEBYE_HUCKEL_POTENTIAL
#define MJOLNIR_DEBYE_HUCKEL_POTENTIAL
#include <mjolnir/core/constants.hpp>
#include <mjolnir/core/System.hpp>
#include <vector>
#include <cmath>

namespace mjolnir
{

template<typename traitsT>
class DebyeHuckelPotential
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef real_type parameter_type;
    typedef std::vector<parameter_type> container_type;

    // r_cutoff = cutoff_ratio * debye_length
    constexpr static real_type cutoff_ratio = 5.5;

  public:
    DebyeHuckelPotential() = default;
    DebyeHuckelPotential(const container_type& charges)
        : charges_(charges)
    {}
    DebyeHuckelPotential(container_type&& charges)
        : charges_(std::move(charges))
    {}
    ~DebyeHuckelPotential() = default;

    std::size_t size() const {return charges_.size();}
    void resize(const std::size_t i){charges_.resize(i);}
    void reserve(const std::size_t i){charges_.reserve(i);}
    void clear(){charges_.clear();}

    real_type&       operator[](std::size_t i)       noexcept {return charges_[i];}
    real_type const& operator[](std::size_t i) const noexcept {return charges_[i];}
    real_type&       at(std::size_t i)       {return charges_.at(i);}
    real_type const& at(std::size_t i) const {return charges_.at(i);}

    void update(const System<traitsT>& sys)
    {
        if(temperature_ == sys.temperature() &&
           ionic_strength_ == sys.ionic_strength()) return;

        temperature_    = sys.temperature();
        ionic_strength_ = sys.ionic_strength();
        calc_parameters();
        return;
    }

    real_type potential(const std::size_t i, const std::size_t j,
                        const real_type r) const;

    real_type derivative(const std::size_t i, const std::size_t j,
                         const real_type r) const;

    real_type max_cutoff_length() const;

    std::string name() const noexcept {return "DebyeHuckel";}

  private:

    void      calc_parameters();
    real_type calc_dielectric_water(const real_type T, const real_type c) const;

  private:

    real_type temperature_;
    real_type ionic_strength_;
    real_type inv_4_pi_eps0_epsk_;
    real_type debye_length_;
    real_type inv_debye_length_;
    container_type charges_;
};


template<typename traitsT>
inline typename DebyeHuckelPotential<traitsT>::real_type
DebyeHuckelPotential<traitsT>::potential(
        const std::size_t i, const std::size_t j, const real_type r) const
{
    return inv_4_pi_eps0_epsk_ * charges_[i] * charges_[j] *
           std::exp(-r * inv_debye_length_);
}

template<typename traitsT>
inline typename DebyeHuckelPotential<traitsT>::real_type
DebyeHuckelPotential<traitsT>::derivative(
        const std::size_t i, const std::size_t j, const real_type r) const
{
    return inv_4_pi_eps0_epsk_ * inv_debye_length_ * charges_[i] * charges_[j] *
           std::exp(-r * inv_debye_length_) / (r * r);
}

template<typename traitsT>
inline typename DebyeHuckelPotential<traitsT>::real_type
DebyeHuckelPotential<traitsT>::max_cutoff_length() const
{
    return debye_length_ * cutoff_ratio;
}

template<typename traitsT>
void DebyeHuckelPotential<traitsT>::calc_parameters()
{
    const real_type kB   = physics<real_type>::kB;
    const real_type e    = physics<real_type>::e;
    const real_type NA   = physics<real_type>::NA;
    const real_type eps0 = physics<real_type>::vacuum_permittivity;
    const real_type epsk = calc_dielectric_water(temperature_, ionic_strength_);
    const real_type pi   = constants<real_type>::pi;

    this->inv_4_pi_eps0_epsk_ = 1. / (4 * pi * eps0 * epsk);
    this->debye_length_ = std::sqrt(eps0 * epsk * kB * temperature_ /
                                    (2 * NA * e * e * ionic_strength_));
    this->inv_debye_length_ = 1. / this->debye_length_;
    return;
}

template<typename traitsT>
typename DebyeHuckelPotential<traitsT>::real_type
DebyeHuckelPotential<traitsT>::calc_dielectric_water(
        const real_type T, const real_type c) const
{
    return (249.4 - 0.788 * T + 7.2e-4 * T * T) *
           (1. - 2.551 * c + 5.151e-2 * c * c - 6.889e-3 * c * c * c);
}

} // mjolnir
#endif /* MJOLNIR_DEBYE_HUCKEL_POTENTIAL */
