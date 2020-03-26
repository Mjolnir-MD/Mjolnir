#ifndef MJOLNIR_POTENTIAL_LOCAL_COSINE_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_LOCAL_COSINE_POTENTIAL_HPP
#include <mjolnir/math/math.hpp>
#include <limits>

namespace mjolnir
{
template<typename T> class System;

// Cosine potential for dihedral interaction.
// V(r)  = k * (1 + cos(n(r - r0)))
// dV/dr = -k * n * sin(n(r - r0))
template<typename realT>
class CosinePotential
{
  public:
    using real_type = realT;

  public:
    CosinePotential(
        const real_type k, const std::int32_t n, const real_type v0) noexcept
        : k_(k), n_(n), v0_(v0)
    {}
    ~CosinePotential() = default;

    real_type potential(const real_type v) const noexcept
    {
        const real_type dv = v - v0_;
        return k_ * (1 + std::cos(n_ * dv));
    }

    real_type derivative(const real_type v) const noexcept
    {
        const real_type dv = v - v0_;
        return -k_ * n_ * std::sin(n_ * dv);
    }

    template<typename T>
    void initialize(const System<T>&) const noexcept {return;}

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() noexcept {return "Cosine";}

    real_type    k() const noexcept {return k_;}
    std::int32_t n() const noexcept {return n_;}
    real_type   v0() const noexcept {return v0_;}

    real_type cutoff() const noexcept // no cutoff exists.
    {return std::numeric_limits<real_type>::infinity();}

  private:

    real_type    k_;
    std::int32_t n_;
    real_type   v0_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class CosinePotential<double>;
extern template class CosinePotential<float>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif /* MJOLNIR_POTENTIAL_LOCAL_COSINE_POTENTIAL_HPP */
