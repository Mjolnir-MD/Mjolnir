#ifndef MJOLNIR_POTENTIAL_LOCAL_FLEXIBLE_LOCAL_DIHEDRAL_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_LOCAL_FLEXIBLE_LOCAL_DIHEDRAL_POTENTIAL_HPP
#include <mjolnir/math/constants.hpp>
#include <array>
#include <algorithm>
#include <limits>
#include <cmath>

namespace mjolnir
{

template<typename T> class System;

// Flexible Local Dihedral potential (T. Terakawa and S. Takada Biophys J 2011)
// used in the AICG2+ protein model.
// NOTE: term := {C, k_1^s, k_1^c, k_2^s, k_2^c, k_3^s, k_3^c}
template<typename realT>
class FlexibleLocalDihedralPotential
{
  public:
    using real_type = realT;

  public:

    FlexibleLocalDihedralPotential(const real_type k,
            const std::array<real_type, 7>& term) noexcept
        : k_(k), term_(term)
    {
        real_type phi = -math::constants<real_type>::pi();
        this->min_energy = std::numeric_limits<real_type>::max();
        while(phi < math::constants<real_type>::pi())
        {
            this->min_energy = std::min(this->min_energy,
                    term_[1] * std::cos(  phi) + term_[2] * std::sin(  phi) +
                    term_[3] * std::cos(2*phi) + term_[4] * std::sin(2*phi) +
                    term_[5] * std::cos(3*phi) + term_[6] * std::sin(3*phi) +
                    term_[0]);
            phi += 1e-4;
        }
    }
    ~FlexibleLocalDihedralPotential() = default;

    real_type potential(const real_type phi) const noexcept
    {
        const real_type sin1 = std::sin(phi); // expecting compiler replaces
        const real_type cos1 = std::cos(phi); // this with __sincos.
        const real_type sin_sq = sin1 * sin1;
        const real_type cos_sq = cos1 * cos1;
        const real_type sin2 = 2 * sin1 * cos1;
        const real_type cos2 = cos_sq - sin_sq;
        const real_type sin3 = sin1 * (3 - 4 * sin_sq);
        const real_type cos3 = cos1 * (4 * cos_sq - 3);

        return k_ * (term_[1] * cos1 + term_[2] * sin1 +
                     term_[3] * cos2 + term_[4] * sin2 +
                     term_[5] * cos3 + term_[6] * sin3 +
                     term_[0] - min_energy);
    }

    real_type derivative(const real_type phi) const noexcept
    {
        const real_type sin1 = std::sin(phi); // expecting compiler replaces
        const real_type cos1 = std::cos(phi); // this with __sincos.
        const real_type sin_sq = sin1 * sin1;
        const real_type cos_sq = cos1 * cos1;
        const real_type sin2 = 2 * sin1 * cos1;
        const real_type cos2 = cos_sq - sin_sq;
        const real_type sin3 = sin1 * (3 - 4 * sin_sq);
        const real_type cos3 = cos1 * (4 * cos_sq - 3);

        return k_*(-1*term_[1] * sin1 +   term_[2] * cos1
                   -2*term_[3] * sin2 + 2*term_[4] * cos2
                   -3*term_[5] * sin3 + 3*term_[6] * cos3);
    }

    template<typename T>
    void initialize(const System<T>&) const noexcept {return;}

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() noexcept {return "FlexibleLocalDihedral";}

    real_type                       k()    const noexcept {return k_;}
    std::array<real_type, 7> const& coef() const noexcept {return term_;}

    real_type cutoff() const noexcept // no cutoff exists.
    {return std::numeric_limits<real_type>::infinity();}

  private:
    real_type min_energy;
    real_type k_;
    std::array<real_type, 7> term_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class FlexibleLocalDihedralPotential<double>;
extern template class FlexibleLocalDihedralPotential<float>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif // MJOLNIR_FLEXIBLE_LOCAL_DIHEDRAL_POTENTIAL
