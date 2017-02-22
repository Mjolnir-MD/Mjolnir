#ifndef MJOLNIR_META
#define MJOLNIR_META
#include <type_traits>

namespace mjolnir
{

template<typename ...>
struct _or_;
template<>
struct _or_<> : public std::false_type {};
template<typename T>
struct _or_<T> : public T {};
template<typename T1, typename T2>
struct _or_<T1, T2> : public std::conditional<T1::value, T1, T2>::type {};
template<typename T1, typename T2, typename T3, typename ... Ts>
struct _or_<T1, T2, T3, Ts...> 
: public std::conditional<T1::value, T1, _or_<T2, T3, Ts...>>::type {};

template<typename ...>
struct _and_;
template<>
struct _and_<> : public std::true_type {};
template<typename T>
struct _and_<T> : public T {};
template<typename T1, typename T2>
struct _and_<T1, T2> : public std::conditional<T1::value, T2, T1>::type {};
template<typename T1, typename T2, typename T3, typename ... Ts>
struct _and_<T1, T2, T3, Ts...> 
: public std::conditional<T1::value, _and_<T2, T3, Ts...>, T1>::type {};

template<typename T>
struct _not_ : public std::integral_constant<bool, !T::value>{};

} // mjolnir
#endif /* MJOLNIR_META */
