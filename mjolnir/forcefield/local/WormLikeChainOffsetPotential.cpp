#include <mjolnir/forcefield/local/WormLikeChainOffsetPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
    template class WormLikeChainOffsetPotential<double>;
    template class WormLikeChainOffsetPotential<float>;
} // mjolnir
