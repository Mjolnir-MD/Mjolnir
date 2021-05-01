#include <mjolnir/forcefield/local/MBasinRepulsivePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class MBasinRepulsivePotential<double>;
template class MBasinRepulsivePotential<float>;
} // mjolnir
