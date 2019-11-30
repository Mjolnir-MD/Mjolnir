#include <mjolnir/potential/local/GoContactAttractivePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GoContactAttractivePotential<double>;
template class GoContactAttractivePotential<float>;
} // mjolnir
