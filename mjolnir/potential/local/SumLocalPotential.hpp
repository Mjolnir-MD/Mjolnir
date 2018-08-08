#ifndef MJOLNIR_SUM_LOCAL_POTENTIAL
#define MJOLNIR_SUM_LOCAL_POTENTIAL
#include <tuple>

namespace mjolnir
{
template<typename T> class System;

namespace detail
{

template<std::size_t Idx, std::size_t Last,
         typename realT, typename ... Potentials>
struct sumup_potential_impl
{
    using real_type = realT;

    inline static real_type
    invoke(const real_type x, const std::tuple<Potentials...>& pots) noexcept
    {
        return std::get<Idx>(pots).potential(x) +
            sumup_potential_impl<Idx+1, Last, realT, Potentials...>::invoke(x, pots);
    }
};
template<std::size_t Last, typename realT, typename ... Potentials>
struct sumup_potential_impl<Last, Last, realT, Potentials ...>
{
    using real_type = realT;
    inline static real_type
    invoke(const real_type x, const std::tuple<Potentials...>& pots) noexcept
    {
        return 0;
    }
};

template<std::size_t Idx, std::size_t Last,
         typename realT, typename ... Potentials>
struct sumup_derivative_impl
{
    using real_type = realT;
    inline static real_type
    invoke(const real_type x, const std::tuple<Potentials...>& pots) noexcept
    {
        return std::get<Idx>(pots).derivative(x) +
            sumup_derivative_impl<Idx+1, Last, realT,  Potentials...>::invoke(x, pots);
    }
};
template<std::size_t Last, typename realT, typename ... Potentials>
struct sumup_derivative_impl<Last, Last, realT, Potentials ...>
{
    using real_type = realT;
    inline static real_type
    invoke(const real_type x, const std::tuple<Potentials...>& pots) noexcept
    {
        return 0.0;
    }
};

template<std::size_t Idx, std::size_t Last, typename ... Potentials>
struct update_all_impl
{
    template<typename T>
    inline static void
    invoke(const System<T>& sys, const std::tuple<Potentials...>& pots) noexcept
    {
        std::get<Idx>(pots).update(sys);
        return update_all_impl<Idx+1, Last, Potentials...>::invoke(
                sys, pots);
    }
};
template<std::size_t Last, typename ... Potentials>
struct update_all_impl<Last, Last, Potentials...>
{
    template<typename T>
    inline static void
    invoke(const System<T>& sys, const std::tuple<Potentials...>& pots) noexcept
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

template<typename realT, typename ... Potentials>
class SumLocalPotential
{
  public:
    using real_type       = realT;
    using potentials_type = std::tuple<Potentials...>;

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
        return detail::sumup_potential_impl<0, size, real_type, Potentials ...
            >::invoke(x, this->potentials_);
    }

    real_type derivative(const real_type x) const noexcept
    {
        return detail::sumup_derivative_impl<0, size, real_type, Potentials ...
            >::invoke(x, this->potentials_);
    }

    template<typename traitsT>
    void update(const System<traitsT>& sys) const noexcept
    {
        return detail::update_all_impl<0, size, Potentials...>::invoke(
                sys, this->potentials_);
    }

    static std::string name()
    {
        return detail::collect_name_impl<0, size, Potentials...>::invoke();
    }

  private:
    potentials_type potentials_;
};
template<typename realT, typename ... Potentials>
constexpr std::size_t SumLocalPotential<realT, Potentials...>::size;

} // mjolnir
#endif//MJOLNIR_SUM_LOCAL_POTENTIAL
