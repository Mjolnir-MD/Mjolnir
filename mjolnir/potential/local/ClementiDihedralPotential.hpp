#ifndef MJOLNIR_CLEMENTI_DIHEDRAL_POTENTIAL
#define MJOLNIR_CLEMENTI_DIHEDRAL_POTENTIAL
#include <cmath>

namespace mjolnir
{
template<typename T> class System;

// triple-cosine potential used as the dihedral term in Clementi's off-lattice
// Go-like protein model (Clement et al., 2000).
//  V(phi) = k1 * (1-cos(phi-phi0)) + k3 * (1-cos(3(phi-phi0)))
// dV/dphi = k1 *    sin(phi-phi0)  + k3 * 3 * sin(3(phi-phi0))
template<typename realT>
class ClementiDihedralPotential
{
  public:
    using real_type = realT;

  public:
    ClementiDihedralPotential(
        const real_type k1, const real_type k3, const real_type v0)
        : k1_(k1), k3_(k3), v0_(v0)
    {}
    ~ClementiDihedralPotential() = default;

    real_type potential(const real_type val) const noexcept
    {
        const real_type dphi = val - v0_;
        const real_type cos1 = std::cos(dphi);
        const real_type cos3 = cos1 * (real_type(4.0) * cos1 * cos1 - real_type(3.0));
        return this->k1_ * (real_type(1.0) - cos1) + k3_ * (real_type(1.0) - cos3);
    }

    real_type derivative(const real_type val) const noexcept
    {
        const real_type dphi = val - v0_;
        const real_type sin1 = std::sin(dphi);
        const real_type sin3 = sin1 * (real_type(3.0) - real_type(4.0) * sin1 * sin1);

        return this->k1_ * sin1 + real_type(3.0) * this->k3_ * sin3;
    }

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() {return "ClementiDihedral";}

  private:

    real_type k1_;
    real_type k3_;
    real_type v0_;
};

} // mjolnir
#endif /* MJOLNIR_CLEMENTI_DIHEDRAL_POTENTIAL */
