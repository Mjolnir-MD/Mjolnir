#ifndef MJOLNIR_CLEMENTI_DIHEDRAL_POTENTIAL
#define MJOLNIR_CLEMENTI_DIHEDRAL_POTENTIAL
#include <string>
#include <cmath>

namespace mjolnir
{
template<typename T> class System;

/*! @brief Clementi-Go dihedral modified-triple-cosine potential *
 *  V(phi) = k1 * (1-cos(phi-phi0)) + k3 * (1-cos(3(phi-phi0)))  *
 * dV/dphi = k1 * sin(phi-phi0)     + k3 * 3 * sin(3(phi-phi0))  */
template<typename traitsT>
class ClementiDihedralPotential
{
  public:
    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:
    ClementiDihedralPotential(const real_type k1, const real_type k3,
                              const real_type native_val)
        : k1_(k1), k3_(k3), native_val_(native_val)
    {}
    ~ClementiDihedralPotential() = default;

    real_type potential(const real_type val) const noexcept
    {
        const real_type dphi = val - native_val_;
        const real_type cos1 = std::cos(dphi);
        const real_type cos3 = cos1 * (4 * cos1 * cos1 - 3.0);
        return this->k1_ * (1.0 - cos1) + k3_ * (1.0 - cos3);
    }

    real_type derivative(const real_type val) const noexcept
    {
        const real_type dphi = val - native_val_;
        const real_type sin1 = std::sin(dphi);
        const real_type sin3 = sin1 * (3.0 - 4 * sin1 * sin1);

        return this->k1_ * sin1 + 3.0 * this->k3_ * sin3;
    }

    void update(const system_type&, const real_type) const noexcept {return;}

    const char* name() const noexcept {return "ClementiDihedral";}

  private:

    real_type k1_;
    real_type k3_;
    real_type native_val_;
};

} // mjolnir
#endif /* MJOLNIR_CLEMENTI_DIHEDRAL_POTENTIAL */
