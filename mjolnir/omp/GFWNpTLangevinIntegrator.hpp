#ifndef MJOLNIR_OMP_GFW_NPT_LANGEVIN_INTEGRATOR_HPP
#define MJOLNIR_OMP_GFW_NPT_LANGEVIN_INTEGRATOR_HPP
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/ForceField.hpp>
#include <mjolnir/omp/RandomNumberGenerator.hpp>
#include <mjolnir/core/GFWNpTLangevinIntegrator.hpp>

namespace mjolnir
{

// a specialization of BAOAB Langevin integrator for OpenMP implementation
template<typename realT>
class GFWNpTLangevinIntegrator<OpenMPSimulatorTraits<realT, CuboidalPeriodicBoundary>>
{
  public:
    using traits_type     = OpenMPSimulatorTraits<realT, CuboidalPeriodicBoundary>;
    using boundary_type   = typename traits_type::boundary_type;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using system_type     = System<traits_type>;
    using forcefield_type = ForceField<traits_type>;
    using rng_type        = RandomNumberGenerator<traits_type>;

  public:

    GFWNpTLangevinIntegrator(const real_type dt, const real_type chi,
            const coordinate_type& m_cell,     const coordinate_type& gamma_cell,
            const coordinate_type& v_cell_ini, const std::vector<real_type>& gammas)
        : dt_(dt), halfdt_(dt / 2), temperature_(/* dummy = */ -1), chi_(chi),
          P_ins_(math::make_coordinate<coordinate_type>(0, 0, 0)),
          P_ref_(math::make_coordinate<coordinate_type>(0, 0, 0)),
          m_cell_(m_cell), v_cell_(v_cell_ini), gamma_cell_(gamma_cell),
          gammas_(gammas), exp_gamma_dt_(gammas.size()), noise_coeff_(gammas.size())
    {}
    ~GFWNpTLangevinIntegrator() = default;

    void initialize(system_type& sys, forcefield_type& ff, rng_type&)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // calculate parameters for each particles
        this->update(sys);

        // calculate force
#pragma omp parallel for
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
        }
        ff.calc_force(sys);

        // calculate the current pressure using the force calculated here
        const auto h_cell = sys.boundary().upper_bound() - sys.boundary().lower_bound();
        this->P_ins_ = this->calc_pressure(sys, h_cell);

        sys.attribute("P_instantaneous_x") = X(P_ins_);
        sys.attribute("P_instantaneous_y") = Y(P_ins_);
        sys.attribute("P_instantaneous_z") = Z(P_ins_);
        sys.attribute("volume")            = X(h_cell) * Y(h_cell) * Z(h_cell);

        return;
    }

    real_type step(const real_type t, system_type& sys,
                   forcefield_type& ff, rng_type& rng)
    {
        using math::hadamard_product;
        using math::make_coordinate;
        using math::X; using math::Y; using math::Z;

        // --------------------------------------------------------------------
        // preparation

        const auto cell_origin = sys.boundary().lower_bound();

        auto h_cell = sys.boundary().upper_bound() - cell_origin;
        auto det_h  = X(h_cell) * Y(h_cell) * Z(h_cell);
        auto inv_h  = make_coordinate<coordinate_type>(
                        1 / X(h_cell), 1 / Y(h_cell), 1 / Z(h_cell));
        auto P_diff = P_ins_ - P_ref_;

        // --------------------------------------------------------------------
        // update cell velocity

        v_cell_ += hadamard_product(halfdt_ * rm_cell_,
                   hadamard_product(det_h * inv_h, P_diff) - chi_kBT_ * inv_h);

        // --------------------------------------------------------------------
        // update particle velocities

        {
            const auto v_over_h = hadamard_product(-v_cell_, inv_h);
            const auto coef_S = make_coordinate<coordinate_type>(
                    std::exp(halfdt_ * X(v_over_h)),
                    std::exp(halfdt_ * Y(v_over_h)),
                    std::exp(halfdt_ * Z(v_over_h)));

#pragma omp parallel for
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                sys.velocity(i) = this->S(sys.velocity(i),
                                          sys.force(i) * sys.rmass(i),
                                          v_over_h, coef_S);
                sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
            }
        }

        // --------------------------------------------------------------------
        // update cell size, volume and inverse cell matrix

        h_cell  += halfdt_ * v_cell_;
        det_h    = X(h_cell) * Y(h_cell) * Z(h_cell);
        inv_h    = make_coordinate<coordinate_type>(
                    1 / X(h_cell), 1 / Y(h_cell), 1 / Z(h_cell));

        sys.boundary().set_boundary(cell_origin, cell_origin + h_cell);

        // --------------------------------------------------------------------
        // update particle positions

        real_type max_displacement_sq_1 = 0;
        {
            const auto v_over_h = hadamard_product(v_cell_, inv_h);
            const auto coef_S = make_coordinate<coordinate_type>(
                    std::exp(halfdt_ * X(v_over_h)),
                    std::exp(halfdt_ * Y(v_over_h)),
                    std::exp(halfdt_ * Z(v_over_h)));

#pragma omp parallel for reduction(max:max_displacement_sq_1)
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                const auto new_pos = this->S(sys.position(i), sys.velocity(i),
                                             v_over_h, coef_S);

                max_displacement_sq_1 = std::max(max_displacement_sq_1,
                        math::length_sq(sys.position(i) - new_pos));

                sys.position(i) = sys.boundary().adjust_position(new_pos);
            }
        }

        // --------------------------------------------------------------------
        // update cell velocity (ornstein-Uhlenbeck)

        const auto R_cell = this->gen_R(rng);

        X(v_cell_) *= X(exp_gamma_dt_cell_);
        Y(v_cell_) *= Y(exp_gamma_dt_cell_);
        Z(v_cell_) *= Z(exp_gamma_dt_cell_);

        X(v_cell_) += X(noise_coeff_cell_) * X(R_cell);
        Y(v_cell_) += Y(noise_coeff_cell_) * Y(R_cell);
        Z(v_cell_) += Z(noise_coeff_cell_) * Z(R_cell);

        // --------------------------------------------------------------------
        // update particle velocities (Ornstein-Uhlenbeck)

#pragma omp parallel for
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto R     = this->gen_R(rng);
            sys.velocity(i) *= this->exp_gamma_dt_[i];
            sys.velocity(i) += this->noise_coeff_[i] * R;
        }

        // --------------------------------------------------------------------
        // update particle positions


        real_type max_displacement_sq_2 = 0;
        {
            const auto v_over_h = hadamard_product(v_cell_, inv_h);
            const auto coef_S = make_coordinate<coordinate_type>(
                    std::exp(halfdt_ * X(v_over_h)),
                    std::exp(halfdt_ * Y(v_over_h)),
                    std::exp(halfdt_ * Z(v_over_h)));

#pragma omp parallel for reduction(max:max_displacement_sq_2)
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                const auto new_pos = this->S(sys.position(i), sys.velocity(i),
                                             v_over_h, coef_S);

                max_displacement_sq_2 = std::max(max_displacement_sq_2,
                        math::length_sq(sys.position(i) - new_pos));

                sys.position(i) = sys.boundary().adjust_position(new_pos);
            }
        }

        // --------------------------------------------------------------------
        // update cell size

        h_cell += halfdt_ * v_cell_;
        det_h   = X(h_cell) * Y(h_cell) * Z(h_cell);
        inv_h   = make_coordinate<coordinate_type>(
                    1 / X(h_cell), 1 / Y(h_cell), 1 / Z(h_cell));

        sys.boundary().set_boundary(cell_origin, cell_origin + h_cell);

        // --------------------------------------------------------------------
        // calc force

        // update neighbor list; reduce margin, reconstruct the list if needed
        const auto max_displacement = std::sqrt(max_displacement_sq_1) +
                                      std::sqrt(max_displacement_sq_2);
        ff.reduce_margin(2 * max_displacement, sys);

        ff.calc_force(sys);

        // --------------------------------------------------------------------
        // update particle velocities

        {
            const auto v_over_h = hadamard_product(-v_cell_, inv_h);
            const auto coef_S = make_coordinate<coordinate_type>(
                    std::exp(halfdt_ * X(v_over_h)),
                    std::exp(halfdt_ * Y(v_over_h)),
                    std::exp(halfdt_ * Z(v_over_h)));

#pragma omp parallel for
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                sys.velocity(i) = this->S(sys.velocity(i),
                                          sys.force(i) * sys.rmass(i),
                                          v_over_h, coef_S);
            }
        }

        // --------------------------------------------------------------------
        // calc current pressure
        this->P_ins_ = this->calc_pressure(sys, h_cell);

        sys.attribute("P_instantaneous_x") = X(P_ins_);
        sys.attribute("P_instantaneous_y") = Y(P_ins_);
        sys.attribute("P_instantaneous_z") = Z(P_ins_);
        sys.attribute("volume")            = det_h;

        P_diff = P_ins_ - P_ref_;

        // --------------------------------------------------------------------
        // update cell velocities

        v_cell_ += hadamard_product(halfdt_ * rm_cell_,
                   hadamard_product(det_h * inv_h, P_diff) - chi_kBT_ * inv_h);

        return t + dt_;
    }

    void update(const system_type& sys)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // update reference temperature and reference pressure
        this->temperature_    = check_attribute_exists(sys, "temperature");
        math::X(this->P_ref_) = check_attribute_exists(sys, "pressure_x");
        math::Y(this->P_ref_) = check_attribute_exists(sys, "pressure_y");
        math::Z(this->P_ref_) = check_attribute_exists(sys, "pressure_z");

        MJOLNIR_LOG_NOTICE("system temperature is ", this->temperature_);
        MJOLNIR_LOG_NOTICE("system pressure is ", this->P_ref_);

        this->reset_parameters(sys);
        return ;
    }

    real_type                     delta_t()    const noexcept {return dt_;}
    std::vector<real_type> const& parameters() const noexcept {return gammas_;}

  private:

    real_type check_attribute_exists(const system_type& sys, const std::string& attr)
    {
        if(!sys.has_attribute(attr))
        {
            throw_exception<std::out_of_range>("mjolnir::GFWNpTLangevinIntegrator: "
                "It requires reference temperature, but `", attr,
                "` is not found in `system.attribute`.");
        }
        return sys.attribute(attr);
    }

    // reset all the cached parameters.
    void reset_parameters(const system_type& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        using math::X; using math::Y; using math::Z;

        const auto kBT = physics::constants<real_type>::kB() * this->temperature_;

        this->chi_kBT_ = this->chi_ * kBT;
        MJOLNIR_LOG_INFO("chi = ", this->chi_, ", chi_kBT = ", chi_kBT_);

        X(rm_cell_) = real_type(1) / X(m_cell_);
        Y(rm_cell_) = real_type(1) / Y(m_cell_);
        Z(rm_cell_) = real_type(1) / Z(m_cell_);

        // exp(-gamma dt)
        const real_type exp_gamma_dt_x = std::exp(-dt_ * X(gamma_cell_));
        const real_type exp_gamma_dt_y = std::exp(-dt_ * Y(gamma_cell_));
        const real_type exp_gamma_dt_z = std::exp(-dt_ * Z(gamma_cell_));

        // exp(-2gamma dt)
        const real_type exp_2gamma_dt_x = exp_gamma_dt_x * exp_gamma_dt_x;
        const real_type exp_2gamma_dt_y = exp_gamma_dt_y * exp_gamma_dt_y;
        const real_type exp_2gamma_dt_z = exp_gamma_dt_z * exp_gamma_dt_z;

        X(exp_gamma_dt_cell_) = exp_gamma_dt_x;
        Y(exp_gamma_dt_cell_) = exp_gamma_dt_y;
        Z(exp_gamma_dt_cell_) = exp_gamma_dt_z;

        // we use velocity in this class, so devide it by M_cell.
        X(noise_coeff_cell_) = std::sqrt((1 - exp_2gamma_dt_x) * kBT * X(rm_cell_));
        Y(noise_coeff_cell_) = std::sqrt((1 - exp_2gamma_dt_y) * kBT * Y(rm_cell_));
        Z(noise_coeff_cell_) = std::sqrt((1 - exp_2gamma_dt_z) * kBT * Z(rm_cell_));

#pragma omp parallel for
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto gamma    = this->gammas_.at(i);
            const auto gamma_dt = -1 * gamma * this->dt_;
            this->exp_gamma_dt_.at(i) = std::exp(gamma_dt);
            this->noise_coeff_ .at(i) = std::sqrt(
                    kBT * (1 - std::exp(2 * gamma_dt)) * sys.rmass(i));
        }
        return;
    }

    // It calculates only diagonal terms. off-diagnoal terms are not considered.
    coordinate_type calc_pressure(const system_type&    sys,
                                  const coordinate_type h_cell) const noexcept
    {
        const auto det_h  = math::X(h_cell) * math::Y(h_cell) * math::Z(h_cell);
        const auto rdet_h = real_type(1) / det_h;

        real_type Px = 0;
        real_type Py = 0;
        real_type Pz = 0;
#pragma omp parallel for reduction(+:Px) reduction(+:Py) reduction(+:Pz)
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto m = sys.mass(i);
            const auto r = sys.position(i);
            const auto v = sys.velocity(i);
            const auto f = sys.force(i);

            // here the original form is a tensor product. but the diagonal term
            // of the tensor product (cij = ai * bj) becomes hadamard product
            // (cii = ai * bi).

            const auto Pi = m * math::hadamard_product(v, v) +
                            math::hadamard_product(f, r);
            Px += math::X(Pi);
            Py += math::Y(Pi);
            Pz += math::Z(Pi);
        }
        return math::make_coordinate<coordinate_type>(
                Px * rdet_h, Py * rdet_h, Pz * rdet_h);
    }

    // When the off-diagonal terms vanish, the Su and the Sl functions become
    // the same with each other and are drastically simplified.
    // Here it implements diagonal version of Su and Sl.
    coordinate_type S(coordinate_type x0,      const coordinate_type b,
                      const coordinate_type A, const coordinate_type coef_S) const
    {
        using math::X; using math::Y; using math::Z;

        X(x0) *= X(coef_S);
        Y(x0) *= Y(coef_S);
        Z(x0) *= Z(coef_S);

        X(x0) -= (real_type(1) - X(coef_S)) * X(b) / X(A);
        Y(x0) -= (real_type(1) - Y(coef_S)) * Y(b) / Y(A);
        Z(x0) -= (real_type(1) - Z(coef_S)) * Z(b) / Z(A);

        return x0;
    }

    coordinate_type gen_R(rng_type& rng) noexcept
    {
        const auto x = rng.gaussian();
        const auto y = rng.gaussian();
        const auto z = rng.gaussian();
        return math::make_coordinate<coordinate_type>(x, y, z);
    }

  private:

    real_type dt_;
    real_type halfdt_;
    real_type temperature_;
    real_type chi_;
    real_type chi_kBT_;

    // diagonal terms only
    coordinate_type P_ins_;
    coordinate_type P_ref_;
    coordinate_type m_cell_;
    coordinate_type rm_cell_;
    coordinate_type v_cell_;
    coordinate_type gamma_cell_;
    coordinate_type exp_gamma_dt_cell_;
    coordinate_type noise_coeff_cell_;

    std::vector<real_type> gammas_;
    std::vector<real_type> exp_gamma_dt_;
    std::vector<real_type> noise_coeff_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class GFWNpTLangevinIntegrator<OpenMPSimulatorTraits<double, UnlimitedBoundary>>;
extern template class GFWNpTLangevinIntegrator<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>;
extern template class GFWNpTLangevinIntegrator<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class GFWNpTLangevinIntegrator<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif /* MJOLNIR_GFW_NPT_LANGEVIN_INTEGRATOR */
