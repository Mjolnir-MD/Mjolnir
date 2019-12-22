#include <mjolnir/potential/global/InversePowerPotential.hpp>

#ifdef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
    template class InversePowerPotential<double>;
    template class InversePowerPotential<float>;
}// mjolnir
