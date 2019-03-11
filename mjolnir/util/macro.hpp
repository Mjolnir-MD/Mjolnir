#ifndef MJOLNIR_UTIL_MACRO_HPP
#define MJOLNIR_UTIL_MACRO_HPP

// MJOLNIR_STRINGIZE(hoge) -> "hoge"
#define MJOLNIR_STRINGIZE_AUX(x) #x
#define MJOLNIR_STRINGIZE(x)     MJOLNIR_STRINGIZE_AUX(x)

// MJOLNIR_FUNC_NAME, expanded into function name.
//
// __PRETTY_FUNCTION__ is too informative. since it shows the template argument,
// the function name become too long.

#if defined(MJOLNIR_DEBUG)
// debug mode. print more verbose function name.
#  if defined(__GNUC__)
#    define MJOLNIR_FUNC_NAME __PRETTY_FUNCTION__
#  elif defined(_MSC_VER)
#    define MJOLNIR_FUNC_NAME __FUNCSIG__
#  else
#    define MJOLNIR_FUNC_NAME __func__
#  endif
#else // release mode. simplify function names.
#  define MJOLNIR_FUNC_NAME __func__
#endif
// #endif

#endif // MJOLNIR_UTIL_MACRO_HPP
