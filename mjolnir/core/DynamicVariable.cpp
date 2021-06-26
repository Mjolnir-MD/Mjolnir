#include <mjolnir/core/DynamicVariable.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template class DynamicVariableBase<double>;
template class DynamicVariableBase<float >;
template class DynamicVariable<double>;
template class DynamicVariable<float >;
template class DefaultDynamicVariable<double>;
template class DefaultDynamicVariable<float >;
template class PeriodicDynamicVariable<double>;
template class PeriodicDynamicVariable<float >;
template class RepulsiveDynamicVariable<double>;
template class RepulsiveDynamicVariable<float >;
} // mjolnir
