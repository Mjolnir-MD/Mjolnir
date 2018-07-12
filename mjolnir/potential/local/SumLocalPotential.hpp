#ifndef MJOLNIR_SUM_LOCAL_POTENTIAL
#define MJOLNIR_SUM_LOCAL_POTENTIAL
#include <tuple>

namespace mjolnir
{
template<typename T> class System;

namespace detail
{

template<std::size_t Idx, std::size_t Last, typename ... Potentials>
struct sumup_potential_impl
{
    inline static real_type
    invoke(const real_type x, const std::tuple<Potentials...>& pots) noexcept
    {
        return std::get<Idx>(pots).potential(x) +
            sumup_potential_impl<Idx+1, Last, Potentials...>::invoke(x, pots);
    }
}
template<std::size_t Last, typename ... Potentials>
struct sumup_potential_impl<Last, Last, Potentials ...>
{
    inline static real_type
    invoke(const real_type x, const std::tuple<Potentials...>& pots) noexcept
    {
        return 0.0;
    }
}

template<std::size_t Idx, std::size_t Last, typename ... Potentials>
struct sumup_derivative_impl
{
    inline static real_type
    invoke(const real_type x, const std::tuple<Potentials...>& pots) noexcept
    {
        return std::get<Idx>(pots).derivative(x) +
            sumup_derivative_impl<Idx+1, Last, Potentials...>::invoke(x, pots);
    }
}
template<std::size_t Last, typename ... Potentials>
struct sumup_derivative_impl<Last, Last, Potentials ...>
{
    inline static real_type
    invoke(const real_type x, const std::tuple<Potentials...>& pots) noexcept
    {
        return 0.0;
    }
}

template<std::size_t Idx, std::size_t Last,
         typename traitsT, typename ... Potentials>
struct update_all_impl
{
    inline static void
    invoke(const System<traitsT>& sys, const real_type dt,
           const std::tuple<Potentials...>& pots) noexcept
    {
        std::get<Idx>(pots).update(sys, dt);
        return update_all_impl<Idx+1, Last, traitsT, Potentials...
            >::invoke(sys, dt, pots);
    }
}
template<std::size_t Last, typename traitsT, typename ... Potentials>
struct update_all_impl<Last, Last, traitsT, Potentials...>
{
    inline static void
    invoke(const System<traitsT>& sys, const real_type dt,
           const std::tuple<Potentials...>& pots) noexcept
    {
        return;
    }
}

template<std::size_t Idx, std::size_t Last, typename ... Potentials>
struct collect_name_impl
{
    inline static std::string invoke(const std::tuple<Potentials...>& pots)
    {
        return std::string(std::get<Idx>(pots).name()) + std::string("+") +
            collect_name_impl<Idx+1, Last, Potentials...>::invoke(pots);
    }
}
template<std::size_t Last, typename ... Potentials>
struct collect_name_impl<Last, Last, Potentials...>
{
    inline static std::string invoke(const std::tuple<Potentials...>& pots)
    {
        return "";
    }
}
} // detail

template<typename traitsT, typename ... Potentials>
class SumLocalPotential
{
  public:
    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type real_type;
    typedef std::tuple<Potentials...> potentials_type;

    static constexpr std::size_t size = sizeof...(Potentials);

  public:

    SumLocalPotential(const Potentials& ... args)
        : potentials_(args...)
    {}
    SumLocalPotential(Potentials&& ... args)
        : potentials_(std::move(args)...)
    {}
    ~SumLocalPotential() = default;

    real_type potential(const real_type x) const noexcept
    {
        return detail::sumup_potential_impl<0, size, Potentials ...>::invoke(
                0.0, x, this->potentials_);
    }

    real_type derivative(const real_type x) const noexcept
    {
        return detail::sumup_derivative_impl<0, size, Potentials ...>::invoke(
                0.0, x, this->potentials_);
    }

    void update(const system_type& sys, const real_type dt) const noexcept
    {
        return detail::update_all_impl<0, size, traitsT, Potentials...>::invoke(
                sys, dt, this->potentials);
    }

    std::string name() const
    {
        return collect_name_impl<0, size, Potentials...>::invoke(
                this->potentials_);
    }

  private:
    potentials_type potentials_;
};

template<typename traitsT, template<typename>class ... Potentials>
constexpr std::size_t SumLocalPotentials<traitsT, Potentials...>::size;

} // mjolnir
#endif//MJOLNIR_SUM_LOCAL_POTENTIAL
