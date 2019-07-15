#include <mjolnir/potential/global/IgnoreGroup.hpp>
#include <mjolnir/core/Topology.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class IgnoreGroup<Topology::group_id_type>;
} // mjolnir
