#ifndef MJOLNIR_DEBYE_HUCKEL_POTENTIAL
#define MJOLNIR_DEBYE_HUCKEL_POTENTIAL
#include <mjolnir/core/constants.hpp>
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
    DebyeHuckelPotential(const real_type T, const real_type ionic_strength,
                         const container_type& charges)
        : temperature_(T), ionic_strength_(ionic_strength), charges_(charges)
    {
        calc_parameters();
    }
    ~DebyeHuckelPotential() = default;

    void set_charges(const container_type& charges);
    void set_charges(container_type&& charges);

    void set_temperature(const real_type T);
    void set_ionic_strength(const real_type I);

    std::size_t size() const {return charges_.size();}
    void resize(const std::size_t i){charges_.resize(i);}
    void reserve(const std::size_t i){charges_.reserve(i);}
    void clear(){charges_.clear();}

    real_type&       operator[](const std::size_t i)       {return charges_[i];}
    real_type const& operator[](const std::size_t i) const {return charges_[i];}
    real_type&       at(const std::size_t i)       {return charges_.at(i);}
    real_type const& at(const std::size_t i) const {return charges_.at(i);}

    real_type potential(const std::size_t i, const std::size_t j,
                        const real_type r) const;

    real_type derivative(const std::size_t i, const std::size_t j,
                         const real_type r) const;

    real_type max_cutoff_length() const;

    void reset_parameter(const std::string&, const real_type);

  private:

    void calc_parameters();
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
typename DebyeHuckelPotential<traitsT>::real_type
DebyeHuckelPotential<traitsT>::potential(
        const std::size_t i, const std::size_t j, const real_type r) const
{
    return inv_4_pi_eps0_epsk_ * charges_[i] * charges_[j] *
           std::exp(-r * inv_debye_length_);
}

template<typename traitsT>
typename DebyeHuckelPotential<traitsT>::real_type
DebyeHuckelPotential<traitsT>::derivative(
        const std::size_t i, const std::size_t j, const real_type r) const
{
    return inv_4_pi_eps0_epsk_ * inv_debye_length_ * charges_[i] * charges_[j] *
           std::exp(-r * inv_debye_length_) / (r * r);
}

template<typename traitsT>
typename DebyeHuckelPotential<traitsT>::real_type
DebyeHuckelPotential<traitsT>::max_cutoff_length() const
{
    return debye_length_ * cutoff_ratio;
}

template<typename traitsT>
void DebyeHuckelPotential<traitsT>::set_charges(const container_type& charges)
{
    this->charges_ = charges;
    return;
}

template<typename traitsT>
void DebyeHuckelPotential<traitsT>::set_charges(container_type&& charges)
{
    this->charges_ = std::forward<container_type>(charges);
    return;
}

template<typename traitsT>
void DebyeHuckelPotential<traitsT>::set_temperature(const real_type T)
{
    this->temperature_ = T;
    this->calc_parameters();
}

template<typename traitsT>
void DebyeHuckelPotential<traitsT>::set_ionic_strength(const real_type I)
{
    this->ionic_strength_ = I;
    this->calc_parameters();
}

template<typename traitsT>
void DebyeHuckelPotential<traitsT>::calc_parameters()
{
    const real_type kB   = physics<traitsT>::kB;
    const real_type e    = physics<traitsT>::e;
    const real_type NA   = physics<traitsT>::NA;
    const real_type eps0 = physics<traitsT>::vacuum_permittivity;
    const real_type epsk = calc_dielectric_water(temperature_, ionic_strength_);
    const real_type pi   = constants<traitsT>::pi;

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

template<typename traitsT>
void DebyeHuckelPotential<traitsT>::reset_parameter(
        const std::string& name, const real_type val)
{
    if(name == "temperature")
        this->set_temperature(val);
    else if(name == "ionic_strength")
        this->set_ionic_strength(val);
    return;
}


} // mjolnir
#endif /* MJOLNIR_DEBYE_HUCKEL_POTENTIAL */
