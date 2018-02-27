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

template<std::size_t N>
struct index_sequence_generator
{
    typedef typename index_sequence_concatenator<
            typename index_sequence_generator<N-1>::type,
            index_sequence<N>
        >::type type;
};

template<>
struct index_sequence_generator<0>
{
    typedef index_sequence<0> type;
};
} // detail

template<std::size_t N>
using make_index_sequence =
    typename detail::index_sequence_generator<N-1>::type;

} // mjolnir
#endif// MJOLNIR_INDEX_SEQUENCE_H
