#ifndef MJOLNIR_TEST_UTIL_CLEAR_FORCE_HPP
#define MJOLNIR_TEST_UTIL_CLEAR_FORCE_HPP
#include <mjolnir/math/Vector.hpp>
#include <mjolnir/core/System.hpp>

namespace mjolnir
{
namespace test
{

template<typename traitsT>
void clear_force(System<traitsT>& sys)
{
    using matrix33_type = typename traitsT::matrix33_type;

    for(std::size_t idx=0; idx<sys.size(); ++idx)
    {
        math::X(sys.force(idx)) = 0;
        math::Y(sys.force(idx)) = 0;
        math::Z(sys.force(idx)) = 0;
    }
    sys.virial() = matrix33_type(0,0,0, 0,0,0, 0,0,0);

    return;
}

template<typename traitsT>
void clear_everything(System<traitsT>& sys)
{
    using real_type       = typename traitsT::real_type;
    using coordinate_type = typename traitsT::coordinate_type;
    using matrix33_type   = typename traitsT::matrix33_type;

    for(std::size_t idx=0; idx<sys.size(); ++idx)
    {
        sys.mass(idx)  = real_type(1);
        sys.rmass(idx) = real_type(1);
        sys.position(idx) = math::make_coordinate<coordinate_type>(0,0,0);
        sys.velocity(idx) = math::make_coordinate<coordinate_type>(0,0,0);
        sys.force(idx)    = math::make_coordinate<coordinate_type>(0,0,0);
        sys.name(idx)     = "X";
        sys.group(idx)    = "NONE";
    }
    sys.virial() = matrix33_type(0,0,0, 0,0,0, 0,0,0);

    return;
}

} // test
} // mjolnir
#endif// MJOLNIR_TEST_UTIL_CHECK_FORCE_HPP
