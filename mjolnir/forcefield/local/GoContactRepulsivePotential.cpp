#include <mjolnir/forcefield/local/GoContactRepulsivePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class GoContactRepulsivePotential<double>;
template class GoContactRepulsivePotential<float>;
}
