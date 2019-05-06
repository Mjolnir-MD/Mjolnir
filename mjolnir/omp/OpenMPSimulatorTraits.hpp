#ifndef MJOLNIR_OMP_OPENMP_SIMULATOR_TRAITS
#define MJOLNIR_OMP_OPENMP_SIMULATOR_TRAITS
#include <mjolnir/math/Vector.hpp>
#include <omp.h>

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT>
struct OpenMPSimulatorTraits
{
    using real_type       = realT;
    using coordinate_type = math::Vector<real_type, 3>;

    using matrix33_type = math::Matrix<real_type, 3, 3>;
    using matrix44_type = math::Matrix<real_type, 4, 4>;

    using boundary_type = boundaryT<real_type, coordinate_type>;
};

} // mjolnir
#endif /* MJOLNIR_DEFAULT_TRAITS */
