#ifndef MJOLNIR_UTIL_MAKE_UNIQUE_HPP
#define MJOLNIR_UTIL_MAKE_UNIQUE_HPP
#include <memory>

namespace mjolnir
{

template<typename T, typename ... Ts>
inline std::unique_ptr<T> make_unique(Ts&& ... args)
{
    return std::unique_ptr<T>(new T(std::forward<Ts>(args)...));
}

} // mjolnir
#endif // MJOLNIR_UTIL_MAKE_UNIQUE_HPP
