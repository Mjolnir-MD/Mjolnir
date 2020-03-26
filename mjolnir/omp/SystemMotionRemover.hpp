#ifndef MJOLNIR_OMP_SYSTEM_MOTION_REMOVER_HPP
#define MJOLNIR_OMP_SYSTEM_MOTION_REMOVER_HPP
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/core/SystemMotionRemover.hpp>

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT>
class SystemMotionRemover<OpenMPSimulatorTraits<realT, boundaryT>>
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

    void remove(system_type& sys) const noexcept
    {
        real_type E_kinetic_pre(0);
        if(rescale_)
        {
#pragma omp parallel for reduction(+:E_kinetic_pre)
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                E_kinetic_pre += sys.mass(i) * math::length_sq(sys.velocity(i));
            }
        }
        E_kinetic_pre *= 0.5;

        if(translation_)
        {
            real_type m_tot(0);
            real_type trans_x(0);
            real_type trans_y(0);
            real_type trans_z(0);
            // TODO OpenMP user-defined reduction operator

#pragma omp parallel for reduction(+:E_kinetic_pre) \
                         reduction(+:trans_x) reduction(+:trans_y) reduction(+:trans_z)
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                m_tot += sys.mass(i);
                trans_x += sys.mass(i) * math::X(sys.velocity(i));
                trans_y += sys.mass(i) * math::Y(sys.velocity(i));
                trans_z += sys.mass(i) * math::Z(sys.velocity(i));
            }
            const auto translation = math::make_coordinate<coordinate_type>(
                    trans_x / m_tot, trans_y / m_tot, trans_z / m_tot);

#pragma omp parallel for
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                sys.velocity(i) -= translation;
            }
        }

        if(rotation_)
        {
            real_type m_tot(0);
            real_type com_x(0);
            real_type com_y(0);
            real_type com_z(0);

#pragma omp parallel for reduction(+:m_tot) \
                         reduction(+:com_x) reduction(+:com_y) reduction(+:com_z)
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                m_tot += sys.mass(i);
                com_x += sys.mass(i) * math::X(sys.position(i));
                com_y += sys.mass(i) * math::Y(sys.position(i));
                com_z += sys.mass(i) * math::Z(sys.position(i));
            }

            const auto CoM = math::make_coordinate<coordinate_type>(
                    com_x / m_tot,  com_y / m_tot,  com_z / m_tot);

            // total angular momentum
            real_type Lx(0);
            real_type Ly(0);
            real_type Lz(0);

            // inertia tensor
            real_type I00(0);
            real_type I11(0);
            real_type I22(0);
            real_type I01(0);
            real_type I02(0);
            real_type I12(0);

#pragma omp parallel for reduction(+:Lx)  reduction(+:Ly)  reduction(+:Lz)  \
                         reduction(+:I00) reduction(+:I11) reduction(+:I22) \
                         reduction(+:I01) reduction(+:I02) reduction(+:I12)
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                const auto  m = sys.mass(i);
                const auto  r = sys.position(i) - CoM;
                const auto& v = sys.velocity(i);

                const auto rxv = math::cross_product(r, v);

                Lx += m * math::X(rxv);
                Ly += m * math::Y(rxv);
                Lz += m * math::Z(rxv);

                const real_type Ixx = m * math::X(r) * math::X(r);
                const real_type Iyy = m * math::Y(r) * math::Y(r);
                const real_type Izz = m * math::Z(r) * math::Z(r);
                const real_type Ixy = m * math::X(r) * math::Y(r);
                const real_type Iyz = m * math::Y(r) * math::Z(r);
                const real_type Izx = m * math::Z(r) * math::X(r);

                // diagonal
                I00 += Iyy + Izz;
                I11 += Izz + Ixx;
                I22 += Ixx + Iyy;

                // off-diagonal
                I01 -= Ixy;
                I02 -= Izx;
                I12 -= Iyz;
            }

            // total angular momentum
            const auto L = math::make_coordinate<coordinate_type>(Lx, Ly, Lz);

            // inertia tensor
            const matrix33_type I(I00, I01, I02,
                                  I01, I11, I12,
                                  I02, I12, I22);


            const real_type detI = math::determinant(I);
            if(detI != real_type(0))
            {
                const auto invI = math::inverse(I);
                const auto omega = invI * L; // angular velocity
#pragma omp parallel for
                for(std::size_t i=0; i<sys.size(); ++i)
                {
                    sys.velocity(i) -= math::cross_product(omega, sys.position(i) - CoM);
                }
            }
        }

        if(rescale_)
        {
            real_type E_kinetic_post(0);
#pragma omp parallel for reduction(+:E_kinetic_post)
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                E_kinetic_post += sys.mass(i) * math::length_sq(sys.velocity(i));
            }
            E_kinetic_post *= 0.5;

            // mv^2 / 2, scale for the velocity should be sqrt-ed
            const real_type scale = std::sqrt(E_kinetic_pre / E_kinetic_post);
#pragma omp parallel for
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                sys.velocity(i) *= scale;
            }
        }
        return;
    }

  protected:
    bool translation_;
    bool rotation_;
    bool rescale_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class SystemMotionRemover<OpenMPSimulatorTraits<double, UnlimitedBoundary>       >;
extern template class SystemMotionRemover<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class SystemMotionRemover<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class SystemMotionRemover<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif // SEPARATE_BUILD

} // mjolnir
#endif /* MJOLNIR_MOLECULAR_DYNAMICS_SIMULATOR */
