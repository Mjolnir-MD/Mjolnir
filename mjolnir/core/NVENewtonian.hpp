#ifndef MJOLNIR_NVE_NEWTONIAN_INTEGRATOR
#define MJOLNIR_NVE_NEWTONIAN_INTEGRATOR
#include "Integrator.hpp"
#include "BoundaryCondition.hpp"
#include <mjolnir/util/zip_iterator.hpp>
#include <mjolnir/util/make_zip.hpp>

namespace mjolnir
{

template<typename traitsT, typename boundaryT = UnlimitedBoundary<traitsT>>
class NVENewtonian : public Integrator<traitsT>
{
  public:
    typedef traitsT traits_type;
    typedef boundaryT boundary_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:

    NVENewtonian(const time_type dt, const std::size_t number_of_particles)
        : dt_(dt), halfdt_(dt * 0.5), halfdt2_(dt * dt * 0.5),
          acceleration_(number_of_particles)
    {}
    ~NVENewtonian() override = default;

    void initialize(const ParticleContainer<traitsT>& pcon) override;
    time_type step(const time_type time, ParticleContainer<traitsT>& pcon,
                   ForceField<traitsT>& ff) override;

    time_type& delta_t()       override {return dt_;}
    time_type  delta_t() const override {return dt_;}

  private:

    void adjust_translation(ParticleContainer<traitsT>& pcon) const;
    void adjust_rotation(ParticleContainer<traitsT>& pcon) const;

  private:
    time_type dt_;      //!< dt
    time_type halfdt_;  //!< dt/2
    time_type halfdt2_; //!< dt^2/2
    std::vector<coordinate_type> acceleration_;
};

template<typename traitsT, typename boundaryT>
void NVENewtonian<traitsT, boundaryT>::initialize(
        const ParticleContainer<traitsT>& pcon)
{
    for(auto iter = make_zip(pcon.cbegin(), acceleration_.begin());
            iter != make_zip(pcon.cend(), acceleration_.end()); ++iter)
    {
        *get<1>(iter) = get<0>(iter)->force / get<0>(iter)->mass;
    }

    return;
}

// at the initial step, acceleration_ must be initialized
template<typename traitsT, typename boundaryT>
typename NVENewtonian<traitsT, boundaryT>::time_type
NVENewtonian<traitsT, boundaryT>::step(const time_type time,
        ParticleContainer<traitsT>& pcon, ForceField<traitsT>& ff)
{
    // calc r(t+dt)
    for(auto iter = make_zip(pcon.begin(), acceleration_.cbegin());
            iter != make_zip(pcon.end(), acceleration_.cend()); ++iter)
    {
        get<0>(iter)->position = boundary_type::adjust_absolute(
            get<0>(iter)->position + dt_ * (get<0>(iter)->velocity) +
            halfdt2_ * (*get<1>(iter)));
        get<0>(iter)->velocity += halfdt_ * (*get<1>(iter));
    }

    // calc f(t+dt)
    ff.calc_force(pcon);

    // calc a(t+dt) and v(t+dt)
    for(auto iter = make_zip(pcon.begin(), acceleration_.begin());
            iter != make_zip(pcon.end(), acceleration_.end()); ++iter)
    {
        // consider cash 1/m
        const coordinate_type acc = get<0>(iter)->force / (get<0>(iter)->mass);
        *get<1>(iter) = acc;
        get<0>(iter)->velocity += halfdt_ * acc;
        get<0>(iter)->force = coordinate_type(0., 0., 0.);
    }

    this->adjust_translation(pcon);
    this->adjust_rotation(pcon);

    return time + dt_;
}

template<typename traitsT, typename boundaryT>
void NVENewtonian<traitsT, boundaryT>::adjust_translation(
        ParticleContainer<traitsT>& pcon) const
{
    coordinate_type translation(0, 0, 0);
    for(auto iter = pcon.begin(); iter != pcon.end(); ++iter)
        translation += iter->velocity;

    translation *= (1. / static_cast<real_type>(pcon.size()));

    for(auto iter = pcon.begin(); iter != pcon.end(); ++iter)
        iter->velocity -= translation;
    return;
}

template<typename traitsT, typename boundaryT>
void NVENewtonian<traitsT, boundaryT>::adjust_rotation(
        ParticleContainer<traitsT>& pcon) const
{
    coordinate_type L(0, 0, 0); // total angular momentum
    for(auto iter = pcon.cbegin(); iter != pcon.cend(); ++iter)
        L += iter->mass * cross_product(iter->position, iter->velocity);

    // TODO 1. define matrix_type in traits type
    //      2. SymmetricMatrix is needed
    Matrix<real_type, 3, 3> I(0,0,0,0,0,0,0,0,0);// inertia tensor
    for(auto iter = pcon.cbegin(); iter != pcon.cend(); ++iter)
    {
        const real_type        m = iter->mass;
        const coordinate_type& r = iter->position;

        const real_type ixx = m * r[0] * r[0];
        const real_type iyy = m * r[1] * r[1];
        const real_type izz = m * r[2] * r[2];
        const real_type ixy = m * r[0] * r[1];
        const real_type iyz = m * r[1] * r[2];
        const real_type izx = m * r[2] * r[1];

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

    // TODO more sophisticated way like LU decomposition is needed
    const coordinate_type omega = inverse(I) * L; // angular velocity

    for(auto iter = pcon.begin(); iter != pcon.end(); ++iter)
        iter->velocity -= cross_product(omega, iter->position);

    return;
}


} // mjolnir
#endif /* MJOLNIR_VELOCITY_VERLET_INTEGRATOR */
