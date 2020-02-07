#ifndef MJOLNIR_UTIL_MACRO_HPP
#define MJOLNIR_UTIL_MACRO_HPP

// MJOLNIR_STRINGIZE(hoge) -> "hoge"
#define MJOLNIR_STRINGIZE_AUX(x) #x
#define MJOLNIR_STRINGIZE(x)     MJOLNIR_STRINGIZE_AUX(x)

// MJOLNIR_FUNC_NAME, expanded into function name.
#if defined(__GNUC__)
#  define MJOLNIR_FUNC_NAME __PRETTY_FUNCTION__
#elif defined(_MSC_VER)
#  define MJOLNIR_FUNC_NAME __FUNCSIG__
#else
#  define MJOLNIR_FUNC_NAME __func__
#endif

// check the architecture uses little endian or big endian.
// GCC and clang have __BYTE_ORDER__ macro. The value of the macro is equal to
// __ORDER_BIG_ENDIAN__ if the platform uses big endian, or
// __ORDER_BIG_ENDIAN__ if the platform uses little endian.
//
// MSVC seems not to have this macro. But, practically, we can assume all the
// windows machines use little endian.
//
// For other platform, define MJOLNIR_WITH_LITTLE_ENDIAN or BIG_ENDIAN manually.
// It does not overwrite the status.
#if !defined(MJOLNIR_WITH_LITTLE_ENDIAN) && !defined(MJOLNIR_WITH_BIG_ENDIAN)
#  ifdef __GNUC__
#    if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
#      define MJOLNIR_WITH_BIG_ENDIAN
#    elif __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#      define MJOLNIR_WITH_LITTLE_ENDIAN
#    else
#      error "Mjolnir supports only little or big endian."
#    endif
#  elif defined(_MSC_VER)
#      pragma message("Now compiling on Windows, assuming little endian...")
#      define MJOLNIR_WITH_LITTLE_ENDIAN
#  else
#    error "Unknown platform. Please define MJOLNIR_WITH_LITTLE_ENDIAN or MJOLNIR_WITH_BIG_ENDIAN"
#  endif
#endif

#endif // MJOLNIR_UTIL_MACRO_HPP
