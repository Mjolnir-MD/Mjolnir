#include <mjolnir/potential/local/GoContactPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GoContactPotential<double>;
template class GoContactPotential<float>;
} // mjolnir
