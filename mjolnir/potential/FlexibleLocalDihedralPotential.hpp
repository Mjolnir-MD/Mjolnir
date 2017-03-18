#ifndef MJOLNIR_FLEXIBLE_LOCAL_DIHEDRAL_POTENTIAL
#define MJOLNIR_FLEXIBLE_LOCAL_DIHEDRAL_POTENTIAL
#include <mjolnir/core/constants.hpp>
#include <array>
#include <algorithm>
#include <limits>
#include <cmath>

namespace mjolnir
{

/* @brief Flexible Local dihedral potential
 * V(phi) = \sum_m^3 k_m^s sin(mphi) + \sum_m^3 k_n^c cos(nphi) + C
 * NOTE: term := {C, k_1^s, k_1^c, k_2^s, k_2^c, k_3^s, k_3^c} */
template<typename traitsT>
class FlexibleLocalDihedralPotential
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:

    FlexibleLocalDihedralPotential(const real_type k,
            const std::array<real_type, 7>& term)
        : k_(k), term_(term)
    {
        real_type phi = -constants<traits_type>::pi;
        this->min_energy = std::numeric_limits<real_type>::max();
        while(phi < constants<traitsT>::pi)
        {
            this->min_energy = std::min(this->min_energy,
                    term_[1] * std::cos(  phi) + term_[2] * std::sin(  phi) +
                    term_[3] * std::cos(2*phi) + term_[4] * std::sin(2*phi) +
                    term_[5] * std::cos(3*phi) + term_[6] * std::sin(3*phi) +
                    term_[0]);
            phi += 1e-4;
        }
    }

    real_type potential(const real_type phi) const
    {
        return k_ * (term_[1] * std::cos(  phi) + term_[2] * std::sin(  phi) +
                     term_[3] * std::cos(2*phi) + term_[4] * std::sin(2*phi) +
                     term_[5] * std::cos(3*phi) + term_[6] * std::sin(3*phi) +
                     term_[0] - min_energy);
    }

    real_type derivative(const real_type phi) const
    {
        return k_*(-1*term_[1] * std::sin(  phi) +   term_[2] * std::cos(  phi)
                   -2*term_[3] * std::sin(2*phi) + 2*term_[4] * std::cos(2*phi)
                   -3*term_[5] * std::sin(3*phi) + 3*term_[6] * std::cos(3*phi));
    }

    void reset_parameter(const std::string&, const real_type){return;}

  private:
    real_type min_energy;
    const real_type k_;
    const std::array<real_type, 7> term_;
};


} // mjolnir
#endif // MJOLNIR_FLEXIBLE_LOCAL_DIHEDRAL_POTENTIAL
