#include <mjolnir/potential/global/IgnoreMolecule.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class IgnoreMoleculeBase<Topology::molecule_id_type>;
template class IgnoreNothing     <Topology::molecule_id_type>;
template class IgnoreSelf        <Topology::molecule_id_type>;
template class IgnoreOthers      <Topology::molecule_id_type>;
template class IgnoreMolecule    <Topology::molecule_id_type>;
} // mjolnir
