#ifndef MJOLNIR_IS_ALL
#define MJOLNIR_IS_ALL

namespace mjolnir
{

template<template<typename T1, typename T2> class T_condition,
         typename T_ref, typename T_front, typename ... T_args>
struct is_all
    : std::integral_constant<bool, T_condition<T_ref, T_front>::value &&
                                   is_all<T_condition, T_ref, T_args...>::value>
{};

template<template<typename T1, typename T2> class T_condition,
         typename T_ref, typename T_last>
struct is_all<T_condition, T_ref, T_last>
    : std::integral_constant<bool, T_condition<T_ref, T_last>::value>
{};

} // mjolnir

#endif/* MJOLNIR_IS_ALL */
