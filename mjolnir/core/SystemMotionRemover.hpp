#ifndef MJOLNIR_CORE_SYSTEM_MOTION_REMOVER_HPP
#define MJOLNIR_CORE_SYSTEM_MOTION_REMOVER_HPP
#include <mjolnir/math/math.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

template<typename traitsT>
class SystemMotionRemover
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using matrix33_type   = typename traits_type::matrix33_type;
    using system_type     = System<traits_type>;

    SystemMotionRemover(const bool translation, const bool rotation,
                        const bool rescale)
        : translation_(translation), rotation_(rotation), rescale_(rescale)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        if(!translation && !rotation && rescale)
        {
            MJOLNIR_LOG_ERROR("Invalid input: `translation` and `rotation` are "
                              "not turned on, even though `rescale` is turned on.");
        }
    }

    void remove(system_type& sys) const noexcept;

  protected:
    bool translation_;
    bool rotation_;
    bool rescale_;
};

template<typename traitsT>
inline void SystemMotionRemover<traitsT>::remove(system_type& sys) const noexcept
{
    real_type E_kinetic_pre(0);
    if(rescale_)
    {
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            E_kinetic_pre += sys.mass(i) * math::length_sq(sys.velocity(i));
        }
        E_kinetic_pre *= 0.5;
    }

    if(translation_)
    {
        real_type m_tot(0);
        auto translation = math::make_coordinate<coordinate_type>(0, 0, 0);
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            m_tot += sys.mass(i);
            translation += sys.mass(i) * sys.velocity(i);
        }
        translation *= real_type(1) / m_tot;

        for(std::size_t i=0; i<sys.size(); ++i)
        {
            sys.velocity(i) -= translation;
        }
    }

    if(rotation_)
    {
        real_type m_tot(0);
        auto CoM = math::make_coordinate<coordinate_type>(0, 0, 0);
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            m_tot += sys.mass(i);
            CoM   += sys.mass(i) * sys.position(i);
        }
        CoM *= real_type(1) / m_tot;

        // total angular momentum
        auto L = math::make_coordinate<coordinate_type>(0, 0, 0);

        // inertia tensor
        matrix33_type I(0, 0, 0,
                        0, 0, 0,
                        0, 0, 0);

        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto  m = sys.mass(i);
            const auto  r = sys.position(i) - CoM;
            const auto& v = sys.velocity(i);

            L += m * math::cross_product(r, v);

            const real_type Ixx = m * math::X(r) * math::X(r);
            const real_type Iyy = m * math::Y(r) * math::Y(r);
            const real_type Izz = m * math::Z(r) * math::Z(r);
            const real_type Ixy = m * math::X(r) * math::Y(r);
            const real_type Iyz = m * math::Y(r) * math::Z(r);
            const real_type Izx = m * math::Z(r) * math::X(r);

            // diagonal
            I(0, 0) += Iyy + Izz;
            I(1, 1) += Izz + Ixx;
            I(2, 2) += Ixx + Iyy;

            // off-diagonal
            I(0, 1) -= Ixy;
            I(0, 2) -= Izx;
            I(1, 2) -= Iyz;
        }
        I(1, 0) = I(0, 1);
        I(2, 0) = I(0, 2);
        I(2, 1) = I(1, 2);

        const real_type detI = math::determinant(I);
        if(detI != real_type(0))
        {
            const auto invI = math::inverse(I, detI);
            const auto omega = invI * L; // angular velocity
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                sys.velocity(i) -= math::cross_product(omega, sys.position(i) - CoM);
            }
        }
    }

    if(rescale_)
    {
        real_type E_kinetic_post(0);
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            E_kinetic_post += sys.mass(i) * math::length_sq(sys.velocity(i));
        }
        E_kinetic_post *= 0.5;

        // mv^2 / 2, scale for the velocity should be sqrt-ed
        const real_type scale = std::sqrt(E_kinetic_pre / E_kinetic_post);
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            sys.velocity(i) *= scale;
        }
    }
    return;
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class SystemMotionRemover<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class SystemMotionRemover<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class SystemMotionRemover<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class SystemMotionRemover<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif // SEPARATE_BUILD

} // mjolnir
#endif /* MJOLNIR_MOLECULAR_DYNAMICS_SIMULATOR */
