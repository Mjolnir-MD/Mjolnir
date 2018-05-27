#ifndef MJOLNIR_CORE_SYSTEM_MOTION_REMOVER
#define MJOLNIR_CORE_SYSTEM_MOTION_REMOVER
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/Matrix.hpp>
#include <mjolnir/math/Vector.hpp>
#include <algorithm>

namespace mjolnir
{

template<typename Translation, typename Rotation>
struct SystemMotionRemover;

// remove nothing.
template<>
struct SystemMotionRemover<std::false_type, std::false_type>
{
    template<typename traitsT>
    inline static void invoke(System<traitsT>& sys) noexcept {return;}
};

// remove translation only.
template<>
struct SystemMotionRemover<std::true_type, std::false_type>
{
    template<typename traitsT>
    static void invoke(System<traitsT>& sys) noexcept
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
};

// remove rotation only.
template<>
struct SystemMotionRemover<std::false_type, std::true_type>
{
    template<typename traitsT>
    static void invoke(System<traitsT>& sys) noexcept
    {
        typedef typename traitsT::real_type       real_type;
        typedef typename traitsT::coordinate_type coordinate_type;
        typedef typename traitsT::template matrix_type<3,3> matrix33_type;

        coordinate_type L(0, 0, 0); // total angular momentum
        for(const auto& p : sys)
        {
            L += p.mass * cross_product(p.position, p.velocity);
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

            I(0, 0) +=  iyy + izz;
            I(0, 1) += -ixy;
            I(0, 2) += -izx;
            I(1, 1) +=  izz + ixx;
            I(1, 2) += -iyz;
            I(2, 2) +=  ixx + iyy;
        }
        I(1, 0) = I(0, 1);
        I(2, 0) = I(0, 2);
        I(2, 1) = I(1, 2);

        const real_type det = determinant(I);
        if(det != 0.0)
        {
            const coordinate_type omega = inverse(I, det) * L; // angular velocity
            for(auto& particle : sys)
            {
                particle.velocity -= cross_product(omega, particle.position);
            }
        }
        return;
    }
};

// remove both translation and rotation.
template<>
struct SystemMotionRemover<std::true_type, std::true_type>
{
    template<typename traitsT>
    static void invoke(System<traitsT>& sys)  noexcept
    {
        typedef typename traitsT::real_type       real_type;
        typedef typename traitsT::coordinate_type coordinate_type;
        typedef typename traitsT::template matrix_type<3,3> matrix33_type;

        coordinate_type trans(0, 0, 0); // total translation
        coordinate_type L    (0, 0, 0); // total angular momentum
        for(const auto& p : sys)
        {
            trans += p.velocity;
            L     += p.mass * cross_product(p.position, p.velocity);
        }
        trans *= (1. / static_cast<real_type>(sys.size()));

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

        const real_type det = determinant(I);
        if(det != 0.0)
        {
            const coordinate_type omega = inverse(I, det) * L; // angular velocity
            for(auto& particle : sys)
            {
                particle.velocity -= cross_product(omega, particle.position);
                particle.velocity -= trans;
            }
        }
        else
        {
            for(auto& particle : sys)
            {
                particle.velocity -= trans;
            }
        }
        return;
    }
};

} // mjolnir
#endif// MJOLNIR_CORE_SYSTEM_MOTION_REMOVER
