#ifndef MJOLNIR_SUM_LOCAL_POTENTIAL
#define MJOLNIR_SUM_LOCAL_POTENTIAL
#include <tuple>

namespace mjolnir
{
template<typename T> class System;

namespace detail
{

template<std::size_t Idx, std::size_t Last,
         typename traitsT, typename ... Potentials>
struct sumup_potential_impl
{
    using real_type = typename traitsT::real_type;
    inline static real_type
    invoke(const real_type x, const std::tuple<Potentials...>& pots) noexcept
    {
        return std::get<Idx>(pots).potential(x) +
            sumup_potential_impl<Idx+1, Last, traitsT, Potentials...>::invoke(x, pots);
    }
};
template<std::size_t Last, typename traitsT, typename ... Potentials>
struct sumup_potential_impl<Last, Last, traitsT, Potentials ...>
{
    using real_type = typename traitsT::real_type;
    inline static real_type
    invoke(const real_type x, const std::tuple<Potentials...>& pots) noexcept
    {
        return 0.0;
    }
};

template<std::size_t Idx, std::size_t Last,
         typename traitsT, typename ... Potentials>
struct sumup_derivative_impl
{
    using real_type = typename traitsT::real_type;
    inline static real_type
    invoke(const real_type x, const std::tuple<Potentials...>& pots) noexcept
    {
        return std::get<Idx>(pots).derivative(x) +
            sumup_derivative_impl<Idx+1, Last, traitsT,  Potentials...>::invoke(x, pots);
    }
};
template<std::size_t Last, typename traitsT, typename ... Potentials>
struct sumup_derivative_impl<Last, Last, traitsT, Potentials ...>
{
    using real_type = typename traitsT::real_type;
    inline static real_type
    invoke(const real_type x, const std::tuple<Potentials...>& pots) noexcept
    {
        return 0.0;
    }
};

template<std::size_t Idx, std::size_t Last,
         typename traitsT, typename ... Potentials>
struct update_all_impl
{
    using real_type = typename traitsT::real_type;
    inline static void
    invoke(const System<traitsT>& sys,
           const std::tuple<Potentials...>& pots) noexcept
    {
        std::get<Idx>(pots).update(sys);
        return update_all_impl<Idx+1, Last, traitsT, Potentials...>::invoke(
                sys, pots);
    }
};
template<std::size_t Last, typename traitsT, typename ... Potentials>
struct update_all_impl<Last, Last, traitsT, Potentials...>
{
    using real_type = typename traitsT::real_type;
    inline static void
    invoke(const System<traitsT>& sys,
           const std::tuple<Potentials...>& pots) noexcept
    {
        return;
    }
};

template<std::size_t Idx, std::size_t Last, typename ... Potentials>
struct collect_name_impl
{
    inline static std::string invoke()
    {
        return std::string(std::tuple_element<Idx, std::tuple<Potentials...>
            >::type::name()) + std::string("+") +
            collect_name_impl<Idx+1, Last, Potentials...>::invoke();
    }
};
template<std::size_t Last, typename ... Potentials>
struct collect_name_impl<Last, Last, Potentials...>
{
    inline static std::string invoke()
    {
        return "";
    }
};
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
        return detail::sumup_potential_impl<0, size, traitsT, Potentials ...
            >::invoke(x, this->potentials_);
    }

    real_type derivative(const real_type x) const noexcept
    {
        return detail::sumup_derivative_impl<0, size, traitsT, Potentials ...
            >::invoke(x, this->potentials_);
    }

    void update(const system_type& sys) const noexcept
    {
        return detail::update_all_impl<0, size, traitsT, Potentials...>::invoke(
                sys, this->potentials_);
    }

    static std::string name()
    {
        return detail::collect_name_impl<0, size, Potentials...>::invoke();
    }

  private:
    potentials_type potentials_;
};
template<typename traitsT, typename ... Potentials>
constexpr std::size_t SumLocalPotential<traitsT, Potentials...>::size;

} // mjolnir
#endif//MJOLNIR_SUM_LOCAL_POTENTIAL
