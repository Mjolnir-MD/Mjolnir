#ifndef MJOLNIR_NVE_NEWTONIAN_INTEGRATOR
#define MJOLNIR_NVE_NEWTONIAN_INTEGRATOR
#include "System.hpp"
#include <mjolnir/util/zip_iterator.hpp>
#include <mjolnir/util/make_zip.hpp>

namespace mjolnir
{

template<typename traitsT>
class NewtonianStepper
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::boundary_type   boundary_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef System<traitsT> system_type;
    typedef ForceField<traitsT> forcefield_type;

  public:

    NewtonianStepper(const real_type dt)
        : dt_(dt), halfdt_(dt * 0.5), halfdt2_(dt * dt * 0.5)
    {}
    ~NewtonianStepper() = default;

    void initialize(const system_type& sys);

    real_type step(const real_type time, system_type& sys, forcefield_type& ff);

    real_type delta_t() const {return dt_;}
    void set_delta_t(const real_type dt)
    {
        dt_ = dt; halfdt_ = dt * 0.5; halfdt2_ = halfdt_ * dt_;
    }

  private:

    void adjust_translation(system_type& sys) const;
    void adjust_rotation   (system_type& sys) const;

  private:
    real_type dt_;      //!< dt
    real_type halfdt_;  //!< dt/2
    real_type halfdt2_; //!< dt^2/2
    std::vector<coordinate_type> acceleration_;
};

template<typename traitsT>
void NewtonianStepper<traitsT>::initialize(const system_type& sys)
{
    acceleration_.resize(sys.size());

    for(auto iter = make_zip(sys.cbegin(), acceleration_.begin());
            iter != make_zip(sys.cend(),   acceleration_.end()); ++iter)
    {
        *get<1>(iter) = get<0>(iter)->force / get<0>(iter)->mass;
    }
    return;
}

template<typename traitsT>
typename NewtonianStepper<traitsT>::real_type
NewtonianStepper<traitsT>::step(
        const real_type time, system_type& system, forcefield_type& ff)
{
    // calc r(t+dt)
    for(auto iter = make_zip(system.begin(), acceleration_.cbegin());
            iter != make_zip(system.end(), acceleration_.cend()); ++iter)
    {
        get<0>(iter)->position = system.adjust_position(
            get<0>(iter)->position + dt_ * (get<0>(iter)->velocity) +
            halfdt2_ * (*get<1>(iter)));
        get<0>(iter)->velocity += halfdt_ * (*get<1>(iter));
    }

    // calc f(t+dt)
    ff.calc_force(system);

    // calc a(t+dt) and v(t+dt)
    for(auto iter = make_zip(system.begin(), acceleration_.begin());
            iter != make_zip(system.end(), acceleration_.end()); ++iter)
    {
        // consider cash 1/m
        const coordinate_type acc = get<0>(iter)->force / (get<0>(iter)->mass);
        *get<1>(iter) = acc;
        get<0>(iter)->velocity += halfdt_ * acc;
        get<0>(iter)->force = coordinate_type(0., 0., 0.);
    }

    this->adjust_translation(system);
    this->adjust_rotation(system);

    return time + dt_;
}

template<typename traitsT>
void NewtonianStepper<traitsT>::adjust_translation(system_type& sys) const
{
    coordinate_type translation(0, 0, 0);
    for(auto iter = sys.begin(); iter != sys.end(); ++iter)
        translation += iter->velocity;

    translation *= (1. / static_cast<real_type>(sys.size()));

    for(auto iter = sys.begin(); iter != sys.end(); ++iter)
        iter->velocity -= translation;
    return;
}

template<typename traitsT>
void NewtonianStepper<traitsT>::adjust_rotation(system_type& sys) const
{
    coordinate_type L(0, 0, 0); // total angular momentum
    for(auto iter = sys.cbegin(); iter != sys.cend(); ++iter)
        L += iter->mass * cross_product(iter->position, iter->velocity);

    // TODO 1. define matrix_type in traits type
    //      2. SymmetricMatrix is needed
    Matrix<real_type, 3, 3> I(0,0,0,0,0,0,0,0,0);// inertia tensor
    for(auto iter = sys.cbegin(); iter != sys.cend(); ++iter)
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

    for(auto iter = sys.begin(); iter != sys.end(); ++iter)
        iter->velocity -= cross_product(omega, iter->position);

    return;
}


} // mjolnir
#endif /* MJOLNIR_VELOCITY_VERLET_INTEGRATOR */
