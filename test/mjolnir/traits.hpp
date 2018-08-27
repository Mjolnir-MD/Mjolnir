#ifndef MJOLNIR_TRAITS_FOR_TEST_CODE
#define MJOLNIR_TRAITS_FOR_TEST_CODE
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SystemMotionRemover.hpp>

namespace mjolnir
{
namespace test
{

template<typename realT,
         template<typename, typename> class boundaryT = UnlimitedBoundary>
using traits = mjolnir::SimulatorTraits<realT, boundaryT>;

} // test
} // mjolnir
#endif // MJOLNIR_TRAITS_FOR_TEST_CODE
