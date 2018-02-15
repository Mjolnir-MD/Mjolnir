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

    real_type potential(const std::size_t i, const std::size_t j,
                        const real_type r) const noexcept
    {
        return this->inv_4_pi_eps0_epsk_ *
               this->charges_[i] * this->charges_[j] *
               std::exp(-r * this->inv_debye_length_);
    }

    real_type derivative(const std::size_t i, const std::size_t j,
                         const real_type r) const noexcept
    {
        return this->inv_4_pi_eps0_epsk_ * this->inv_debye_length_ *
               this->charges_[i] * this->charges_[j] *
               std::exp(-r * this->inv_debye_length_) / (r * r);
    }

    real_type max_cutoff_length() const noexcept
    {
        return debye_length_ * cutoff_ratio;
    }

    // for temperature/ionic concentration changes...
    void update(const System<traitsT>& sys) noexcept
    {
        if(temperature_    == sys.temperature() &&
           ionic_strength_ == sys.ionic_strength())
        {
            return;
        }

        // TODO: it can be a shared resource, which is better to manage?
        temperature_    = sys.temperature();
        ionic_strength_ = sys.ionic_strength();
        this->calc_parameters();
        return;
    }

    std::string name() const noexcept {return "DebyeHuckel";}

    // access to the parameters
    std::vector<real_type>&       charges()       noexcept {return charges_;}
    std::vector<real_type> const& charges() const noexcept {return charges_;}

  private:

    void calc_parameters() noexcept
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

    real_type calc_dielectric_water(
            const real_type T, const real_type c) const noexcept
    {
        return (249.4 - 0.788 * T + 7.2e-4 * T * T) *
               (1. - 2.551 * c + 5.151e-2 * c * c - 6.889e-3 * c * c * c);
    }

  private:

    real_type temperature_;
    real_type ionic_strength_;
    real_type inv_4_pi_eps0_epsk_;
    real_type debye_length_;
    real_type inv_debye_length_;
    container_type charges_;
};

template<typename traitsT>
constexpr typename DebyeHuckelPotential<traitsT>::real_type
DebyeHuckelPotential<traitsT>::cutoff_ratio;

} // mjolnir
#endif /* MJOLNIR_DEBYE_HUCKEL_POTENTIAL */
