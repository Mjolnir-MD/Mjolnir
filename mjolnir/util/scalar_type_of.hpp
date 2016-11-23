#ifndef MJOLNIR_SCALAR_TYPE_OF
#define MJOLNIR_SCALAR_TYPE_OF
#include <array>

namespace mjolnir
{

template<typename T>
struct scalar_type_of{};

template<typename T>
struct scalar_type_of<std::array<T>> {typedef T type;};

template<typename T, std::size_t N>
struct scalar_type_of<T[N]> {typedef T type;};

template<typename T>
using scalar_type = typename scalar_type_of<T>::type;

}// mjolnir

#endif /* MJOLNIR_SCALAR_TYPE_OF */
