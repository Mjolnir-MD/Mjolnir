#ifndef MJOLNIR_CORE_SIMULATOR_TRAITS_HPP
#define MJOLNIR_CORE_SIMULATOR_TRAITS_HPP
#include <mjolnir/math/Vector.hpp>

#ifdef MJOLNIR_WITH_OPENMP
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#endif

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT>
struct SimulatorTraits
{
    using real_type       = realT;
    using coordinate_type = math::Vector<real_type, 3>;

    using matrix33_type = math::Matrix<real_type, 3, 3>;
    using matrix44_type = math::Matrix<real_type, 4, 4>;

    using boundary_type = boundaryT<real_type, coordinate_type>;
};

template<typename T>
struct is_simulator_traits : std::false_type{};

template<typename realT, template<typename, typename> class boundaryT>
struct is_simulator_traits<SimulatorTraits<realT, boundaryT>>: std::true_type{};

} // mjolnir
#endif /* MJOLNIR_DEFAULT_TRAITS */
