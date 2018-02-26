#ifndef MJOLNIR_INDEX_SEQUENCE_H
#define MJOLNIR_INDEX_SEQUENCE_H
#include <cstddef>

namespace mjolnir
{

template<std::size_t ... vs>
struct index_sequence{};

namespace detail
{
template<typename T1, typename T2>
struct index_sequence_concatenator;

template<std::size_t ... v1, std::size_t ... v2>
struct index_sequence_concatenator<
    index_sequence<v1...>, index_sequence<v2...>>
{
    typedef index_sequence<v1 ..., v2 ...> type;
};

// [0, 4) -> [0, 2) ++ [2, 4) -> [0, 1) ++ [1, 2) ++ [2, 3) ++ [3, 4)
// [0, 3) -> [0, 1) ++ [1, 3) -> [0, 1) ++ [1, 2) ++ [2, 3)
template<std::size_t first, std::size_t last>
struct index_sequence_generator
{
    constexpr static inline std::size_t
    middle(const std::size_t low, const std::size_t high) noexcept
    {
        return low + (high - low) / 2;
    }

    typedef typename index_sequence_concatenator<
            typename index_sequence_generator<first, middle(first, last)>::type,
            typename index_sequence_generator<middle(first, last), last>::type
        >::type type;
};

template<std::size_t N>
struct index_sequence_generator<N, N+1>
{
    typedef index_sequence<N> type;
};
} // detail

template<std::size_t N>
using make_index_sequence =
    typename detail::index_sequence_generator<0, N>::type;

} // mjolnir
#endif// MJOLNIR_INDEX_SEQUENCE_H
