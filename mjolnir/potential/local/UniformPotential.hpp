#ifndef MJOLNIR_POTENTIAL_UNIFORM_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_UNIFORM_POTENTIAL_HPP
#include <limits>

namespace mjolnir
{
template<typename T> class System;

// Uniform potential is for test potentila.
// V(r) = 1
// dV/dr = 0
template<typename realT>
class UniformPotential
{
  public:
    using real_type = realT;

  public:
    UniformPotential(const real_type k) noexcept
        : k_(k)
    {}
    ~UniformPotential() = default;

    real_type potential(const real_type) const noexcept
    {
        return k_;
    }

    real_type derivative(const real_type) const noexcept
    {
        return 0.0;
    }

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() noexcept {return "Uniform";}

    real_type k()      const noexcept {return k_;}
    real_type cutoff() const noexcept {return std::numeric_limits<real_type>::infinity();}

  private:

    real_type k_;
};

} // mjolnir
#endif /* MJOLNIR_UNIFORM_POTENTIAL */
