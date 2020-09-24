#ifndef MJOLNIR_UTIL_BINARY_IO_HPP
#define MJOLNIR_UTIL_BINARY_IO_HPP
#include <istream>
#include <ostream>
#include <utility>
#include <type_traits>

namespace mjolnir
{
namespace detail
{

template<typename T>
T read_bytes_as(std::istream& os)
{
    T v;
    os.read(reinterpret_cast<char*>(std::addressof(v)), sizeof(T));
    return v;
}

template<typename T>
void write_as_bytes(std::ostream& os, const T& v) noexcept
{
    using Type = typename std::remove_reference<T>::type;
    os.write(reinterpret_cast<const char*>(std::addressof(v)), sizeof(Type));
    return;
}

inline void skip_bytes(std::istream& os, std::size_t N)
{
    os.ignore(N);
    return;
}

} // detail
} // mjolnir
#endif// MJOLNIR_UTIL_BINARY_IO_HPP
