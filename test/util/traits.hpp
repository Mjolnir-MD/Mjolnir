#ifndef MJOLNIR_TEST_UTIL_TRAITS_HPP
#define MJOLNIR_TEST_UTIL_TRAITS_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
namespace test
{

template<typename realT,
         template<typename, typename> class boundaryT = UnlimitedBoundary>
using traits = mjolnir::SimulatorTraits<realT, boundaryT>;

} // test
} // mjolnir
#endif // MJOLNIR_TEST_UTIL_TRAITS_HPP
