#ifndef MJOLNIR_CLEMENTI_DIHEDRAL_POTENTIAL
#define MJOLNIR_CLEMENTI_DIHEDRAL_POTENTIAL
#include <string>
#include <cmath>

namespace mjolnir
{

/*! @brief Clementi-Go dihedral modified-triple-cosine potential *
 *  V(phi) = k1 * (1-cos(phi-phi0)) + k3 * (1-cos(3(phi-phi0)))  *
 * dV/dphi = k1 * sin(phi-phi0)     + k3 * 3 * sin(3(phi-phi0))  */
template<typename traitsT>
class ClementiDihedralPotential
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:
    ClementiDihedralPotential(
            const real_type k1, const real_type k3, const real_type native_val)
        : k1_(k1), k3_(k3), native_val_(native_val)
    {}
    ~ClementiDihedralPotential() = default;

    real_type potential(const real_type val) const
    {
        const real_type dphi = val - native_val_;
        return k1_ * (1. - std::cos(dphi)) + k3_ * (1. - std::cos(3. * dphi));
    }

    real_type derivative(const real_type val) const
    {
        const real_type dphi = val - native_val_;
        return k1_ * std::sin(dphi) + 3. * k3_ * std::sin(3. * dphi);
    }

    std::string name() const noexcept {return "ClementiDihedral";}

  private:

    const real_type k1_;
    const real_type k3_;
    const real_type native_val_;
};

} // mjolnir
#endif /* MJOLNIR_CLEMENTI_DIHEDRAL_POTENTIAL */
