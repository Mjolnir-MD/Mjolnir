#ifndef MJOLNIR_TEST_UTIL_APPLY_RANDOM_PERTURBATION_HPP
#define MJOLNIR_TEST_UTIL_APPLY_RANDOM_PERTURBATION_HPP
#include <mjolnir/math/Vector.hpp>
#include <mjolnir/core/System.hpp>

namespace mjolnir
{
namespace test
{

template<typename traitsT, typename RNG>
void apply_random_perturbation(System<traitsT>& sys, RNG& rng,
                 const typename traitsT::real_type dr)
{
    using real_type = typename traitsT::real_type;
    std::uniform_real_distribution<real_type> uni(-dr, dr);

    for(std::size_t idx=0; idx<sys.size(); ++idx)
    {
        math::X(sys.position(idx)) += uni(rng);
        math::Y(sys.position(idx)) += uni(rng);
        math::Z(sys.position(idx)) += uni(rng);
    }
    return;
}

} // test
} // mjolnir
#endif// MJOLNIR_TEST_UTIL_CHECK_FORCE_HPP
