#ifndef MJOLNIR_CORE_SYSTEM_MOTION_REMOVER
#define MJOLNIR_CORE_SYSTEM_MOTION_REMOVER
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/Matrix.hpp>
#include <mjolnir/math/Vector.hpp>
#include <algorithm>

namespace mjolnir
{

template<typename traitsT>
void remove_translation(System<traitsT>& sys)
{
    typename traitsT::coordinate_type translation(0, 0, 0);
    for(const auto& particle : sys)
    {
        translation += particle.velocity;
    }

    translation *= (1. / static_cast<typename traitsT::real_type>(sys.size()));

    for(auto& particle : sys)
    {
        particle.velocity -= translation;
    }
    return;
}

template<typename traitsT>
void remove_rotation(System<traitsT>& sys)
{
    typedef typename traitsT::real_type       real_type;
    typedef typename traitsT::coordinate_type coordinate_type;
    typedef typename traitsT::template matrix_type<3,3> matrix33_type;

    coordinate_type L(0, 0, 0); // total angular momentum
    for(const auto& particle : sys)
    {
        L += particle.mass * cross_product(particle.position, particle.velocity);
    }

    matrix33_type I(0, 0, 0, // inertia tensor
                    0, 0, 0,
                    0, 0, 0);
    for(const auto& particle : sys)
    {
        const real_type        m = particle.mass;
        const coordinate_type& r = particle.position;

        const real_type ixx = m * r[0] * r[0];
        const real_type iyy = m * r[1] * r[1];
        const real_type izz = m * r[2] * r[2];
        const real_type ixy = m * r[0] * r[1];
        const real_type iyz = m * r[1] * r[2];
        const real_type izx = m * r[2] * r[0];

        I(0, 0) += iyy + izz;
        I(0, 1) += -ixy;
        I(0, 2) += -izx;
        I(1, 1) += izz + ixx;
        I(1, 2) += -iyz;
        I(2, 2) += ixx + iyy;
    }
    I(1, 0) = I(0, 1);
    I(2, 0) = I(0, 2);
    I(2, 1) = I(1, 2);

    // TODO more sophisticated way like LU decomposition is needed?
    const coordinate_type omega = inverse(I) * L; // angular velocity

    for(auto& particle : sys)
    {
        particle.velocity -= cross_product(omega, particle.position);
    }
    return;
}


} // mjolnir
#endif// MJOLNIR_CORE_SYSTEM_MOTION_REMOVER
