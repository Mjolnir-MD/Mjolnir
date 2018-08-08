#ifndef MJOLNIR_COSINE_POTENTIAL
#define MJOLNIR_COSINE_POTENTIAL
#include <cmath>

namespace mjolnir
{

template<typename T> class System;

/*! @brief cosine-potential for 3SPN.1 dihedral potential *
 *  V(phi) = k * (1 - cos(phi - phi0))                    *
 * dV(phi) = k * sin(phi - phi0))                         */
template<typename realT>
class CosinePotential
{
  public:
    using real_type = realT;

  public:
    CosinePotential(const real_type k, const real_type phi0)
        : k_(k), phi0_(phi0)
    {}
    ~CosinePotential() = default;

    real_type potential(const real_type phi) const
    {
        return this->k_ *(1 - std::cos(phi - phi0_));
    }

    real_type derivative(const real_type phi) const
    {
        return this->k_ * std::sin(phi - phi0_);
    }

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() noexcept {return "Cosine";}

  private:

    real_type k_;
    real_type phi0_;
};

} // mjolnir
#endif /* MJOLNIR_COSINE_POTENTIAL */
