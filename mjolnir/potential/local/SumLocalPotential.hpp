#ifndef MJOLNIR_POTENTIAL_LOCAL_SUM_LOCAL_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_LOCAL_SUM_LOCAL_POTENTIAL_HPP
#include <utility>
#include <algorithm>

namespace mjolnir
{
template<typename T> class System;

// It formerly supported an arbitrary number of potentials, but after that
// it turned out that this potential is not so heavily used. Considering the
// maintainability, it is good to remove template tricks and keep it trivial
// by restricting a number of related potentials.
template<typename realT, template<typename T> class Potential1,
                         template<typename T> class Potential2>
class SumLocalPotential
{
  public:
    using real_type       = realT;
    using potential1_type = Potential1<realT>;
    using potential2_type = Potential2<realT>;

  public:

    SumLocalPotential(const Potential1<realT>& p1, const Potential2<realT>& p2)
        : potential1_(p1), potential2_(p2)
    {}
    SumLocalPotential(Potential1<realT>&& p1, Potential2<realT>&& p2)
        : potential1_(std::move(p1)), potential2_(std::move(p2))
    {}
    ~SumLocalPotential() = default;

    real_type potential(const real_type x) const noexcept
    {
        return potential1_.potential(x) + potential2_.potential(x);
    }

    real_type derivative(const real_type x) const noexcept
    {
        return potential1_.derivative(x) + potential2_.derivative(x);
    }

    template<typename traitsT>
    void update(const System<traitsT>& sys) const noexcept
    {
        potential1_.update(sys);
        potential2_.update(sys);
        return;
    }

    static std::string name()
    {
        return potential1_type::name() + std::string("+") +
               potential2_type::name();
    }

    real_type cutoff() const noexcept // no cutoff exists.
    {
        return std::max(potential1_.cutoff() + potential2_.cutoff());
    }

  private:

    potential1_type potential1_;
    potential2_type potential2_;
};

} // mjolnir
#endif//MJOLNIR_SUM_LOCAL_POTENTIAL
