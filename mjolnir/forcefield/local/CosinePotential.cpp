#include <mjolnir/forcefield/local/CosinePotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class CosinePotential<double>;
template class CosinePotential<float>;
} // mjolnir
