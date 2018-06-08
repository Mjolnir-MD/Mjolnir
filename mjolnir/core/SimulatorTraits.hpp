#ifndef MJOLNIR_DEFAULT_TRAITS
#define MJOLNIR_DEFAULT_TRAITS
#include <mjolnir/math/Vector.hpp>

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT>
struct SimulatorTraitsBase
{
    typedef realT                                 real_type;
    typedef Vector<real_type, 3>                  coordinate_type;
    template<std::size_t N, std::size_t M>
    using matrix_type = Matrix<real_type, N, M>;
    typedef coordinate_type                       position_type;
    typedef boundaryT<real_type, coordinate_type> boundary_type;
};

} // mjolnir
#endif /* MJOLNIR_DEFAULT_TRAITS */
