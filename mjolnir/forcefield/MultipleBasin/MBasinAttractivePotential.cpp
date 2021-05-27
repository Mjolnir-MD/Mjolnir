#include <mjolnir/forcefield/MultipleBasin/MBasinAttractivePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class MBasinAttractivePotential<double>;
template class MBasinAttractivePotential<float>;
} // mjolnir
