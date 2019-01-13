#ifndef MJOLNIR_DEFAULT_TRAITS
#define MJOLNIR_DEFAULT_TRAITS
#include <mjolnir/math/Vector.hpp>

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT>
struct SimulatorTraits
{
    using real_type       = realT;
    using coordinate_type = Vector<real_type, 3>;
    template<std::size_t N, std::size_t M>
    using matrix_type = Matrix<real_type, N, M>;

    using boundary_type = boundaryT<real_type, coordinate_type>;
};

} // mjolnir
#endif /* MJOLNIR_DEFAULT_TRAITS */
