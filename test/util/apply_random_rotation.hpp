#ifndef MJOLNIR_TEST_UTIL_APPLY_RANDOM_ROTATION_HPP
#define MJOLNIR_TEST_UTIL_APPLY_RANDOM_ROTATION_HPP
#include <mjolnir/math/math.hpp>
#include <mjolnir/core/System.hpp>

namespace mjolnir
{
namespace test
{

// It rotates the whole system in random direction
template<typename traitsT, typename RNG>
void apply_random_rotation(System<traitsT>& sys, RNG& rng)
{
    using real_type     = typename traitsT::real_type;
    using matrix33_type = typename traitsT::matrix33_type;

    constexpr auto pi = math::constants<real_type>::pi();

    std::uniform_real_distribution<real_type> uni(-pi, pi);

    const auto rot_x = uni(rng);
    const auto rot_y = uni(rng);
    const auto rot_z = uni(rng);

    matrix33_type rotm_x(1.0,             0.0,              0.0,
                         0.0, std::cos(rot_x), -std::sin(rot_x),
                         0.0, std::sin(rot_x),  std::cos(rot_x));

    matrix33_type rotm_y( std::cos(rot_y), 0.0,  std::sin(rot_y),
                                      0.0, 1.0,              0.0,
                         -std::sin(rot_y), 0.0,  std::cos(rot_y));

    matrix33_type rotm_z(std::cos(rot_z), -std::sin(rot_z), 0.0,
                         std::sin(rot_z),  std::cos(rot_z), 0.0,
                                     0.0,              0.0, 1.0);

    const auto rotm = rotm_x * rotm_y * rotm_z;

    for(std::size_t idx=0; idx<sys.size(); ++idx)
    {
        sys.position(idx) = rotm * sys.position(idx);
    }

    return;
}

} // test
} // mjolnir
#endif// MJOLNIR_TEST_UTIL_CHECK_FORCE_HPP
