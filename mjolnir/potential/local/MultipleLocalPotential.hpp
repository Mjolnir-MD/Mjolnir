#ifndef MJOLNIR_MULTIPLE_LOCAL_POTENTIAL
#define MJOLNIR_MULTIPLE_LOCAL_POTENTIAL
#include <tuple>

namespace mjolnir
{
template<typename T> class System;

template<typename traitsT, template<typename>class ... Ts>
class MultipleLocalPotential
{
  public:
    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type real_type;
    constexpr static std::size_t multiplicity = sizeof...(Ts);

  public:

    MultipleLocalPotential(Ts<traits_type> const& ... args)
        : potentials_(args...)
    {}
    MultipleLocalPotential(Ts<traits_type>&& ... args)
        : potentials_(std::move(args)...)
    {}
    ~MultipleLocalPotential() = default;

    real_type potential(const real_type x) const
    {
        return accumurate_potentials<0, multiplicity>::invoke(potentials_, x);
    }

    real_type derivative(const real_type x) const
    {
        return accumurate_derivatives<0, multiplicity>::invoke(potentials_, x);
    }

    void update(const system_type&, const real_type) const noexcept {return;}

  private:

    template<std::size_t I, std::size_t N>
    struct accumurate_potentials
    {
        static constexpr real_type
        invoke(const std::tuple<Ts<traits_type>...>& pots, const real_type x)
        {
            return std::get<I>(pots).potential(x) +
                   accumurate_potentials<I+1, N>::invoke(pots, x);
        }
    };
    template<std::size_t N>
    struct accumurate_potentials<N, N>
    {
        static constexpr real_type
        invoke(const std::tuple<Ts<traits_type>...>& pots, const real_type x)
        {
            return 0.;
        }
    };

    template<std::size_t I, std::size_t N>
    struct accumurate_derivatives
    {
        static constexpr real_type
        invoke(const std::tuple<Ts<traits_type>...>& pots, const real_type x)
        {
            return std::get<I>(pots).derivative(x) +
                   accumurate_derivatives<I+1, N>::invoke(pots, x);
        }
    };
    template<std::size_t N>
    struct accumurate_derivatives<N, N>
    {
        static constexpr real_type
        invoke(const std::tuple<Ts<traits_type>...>& pots, const real_type x)
        {
            return 0.;
        }
    };

  private:
    std::tuple<Ts<traits_type>...> potentials_;
};

}
#endif//MJOLNIR_MULTIPLE_LOCAL_POTENTIAL
