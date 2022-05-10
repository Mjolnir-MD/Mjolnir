#include <mjolnir/forcefield/global/WCAPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class WCAPotential<double>;
template class WCAPotential<float >;
} // mjolnir
