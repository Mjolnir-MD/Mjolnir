#include <mjolnir/potential/global/DebyeHuckelPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class DebyeHuckelPotential<double>;
template class DebyeHuckelPotential<float>;
} // mjolnir
