#include <mjolnir/potential/local/WormLikeChainPotential.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
    template class WormLikeChainPotential<double>;
    template class WormLikeChainPotential<float>;
} // mjolnir
