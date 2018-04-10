#ifndef MJOLNIR_UTIL_TYPE_TRAITS_HPP
#define MJOLNIR_UTIL_TYPE_TRAITS_HPP
#include <type_traits>

namespace mjolnir
{

template<typename ...>
struct conjunction : std::true_type{};
template<typename T>
struct conjunction<T> : T
{
    static_assert(std::is_convertible<decltype(T::value), bool>::value,
                  "conjunction<T> requires T::value is convertible to bool");
};
template<typename T, typename ... Ts>
struct conjunction<T, Ts...> :
    std::conditional<static_cast<bool>(T::value), conjunction<Ts...>, T>::type
{};

template<typename ...>
struct disjunction : std::false_type{};
template<typename T>
struct disjunction<T> : T
{
    static_assert(std::is_convertible<decltype(T::value), bool>::value,
                  "disjunction<T> requires T::value is convertible to bool");
};
template<typename T, typename ... Ts>
struct disjunction<T, Ts...> :
    std::conditional<static_cast<bool>(T::value), T, disjunction<Ts...>>::type
{};

template<typename T>
struct negation : std::integral_constant<bool, !static_cast<bool>(T::value)>{};

template<template<typename> class F, typename ... Ts>
struct is_all : conjunction<F<Ts>...>::type {};

} // mjolnir
#endif// MJOLNIR_UTIL_TYPE_TRAITS_HPP
