#include <mjolnir/core/NeighborList.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class NeighborList<empty_t>;
template class NeighborList<float >;
template class NeighborList<double>;
template class NeighborList<std::pair<float , float >>;
template class NeighborList<std::pair<double, double>>;
} // mjolnir
