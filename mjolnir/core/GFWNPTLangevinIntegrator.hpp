#ifndef MJOLNIR_CORE_GFW_NPT_LANGEVIN_INTEGRATOR_HPP
#define MJOLNIR_CORE_GFW_NPT_LANGEVIN_INTEGRATOR_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceFieldBase.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/core/SystemMotionRemover.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/throw_exception.hpp>

namespace mjolnir
{

// isothermal-isobaric Langevin integrator developed by the following paper
// - Xingyu Gao, Jun Fang, and Han Wang, J. Cmem. Phys. 144, 124113 (2016)
//
// It requires Periodic Boundary Condition. If the UnlimitedBoundary is applied,
// it always throws. PBC implementations are below.
//
// The original algorithm allows cell deformation, but here it inhibits shell
// deformation. The PBC cell will be kept always cuboidal. This condition
// corredponds to the case where we choose Mab -> inf (a != b). Under that
// condition, the off-diagonal terms dissappear.
template<typename traitsT>
class GFWNPTLangevinIntegrator
{
  public:
    using traits_type     = traitsT;
    using boundary_type   = typename traits_type::boundary_type;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using matrix33_type   = typename traits_type::matrix33_type;
    using system_type     = System<traits_type>;
    using forcefield_type = std::unique_ptr<ForceFieldBase<traitsT>>;
    using rng_type        = RandomNumberGenerator<traits_type>;

  public:

    GFWNPTLangevinIntegrator(const real_type dt, const real_type chi,
            const coordinate_type& m_cell,     const coordinate_type& gamma_cell,
            const coordinate_type& v_cell_ini, const std::vector<real_type>& gammas)
        : dt_(dt), halfdt_(dt / 2), temperature_(/* dummy = */ -1), chi_(chi),
          P_ins_(math::make_coordinate<coordinate_type>(0, 0, 0)),
          P_ref_(math::make_coordinate<coordinate_type>(0, 0, 0)),
          m_cell_(m_cell), v_cell_(v_cell_ini), gamma_cell_(gamma_cell),
          gammas_(gammas), exp_gamma_dt_(gammas.size()), noise_coeff_(gammas.size())
    {
        if(!is_cuboidal_periodic_boundary<boundary_type>::value)
        {
            throw_exception<std::runtime_error>("GFWNPTLangevinIntegrator: "
                    "periodic boundary condition is required");
        }
    }
    ~GFWNPTLangevinIntegrator() = default;

    void initialize(system_type& sys, forcefield_type& ff, rng_type&)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        if(!ff->constraint().empty())
        {
            MJOLNIR_LOG_WARN("GFW NPT langevin integrator does not support constraint"
                " forcefield. [[forcefields.constraint]] will be ignored.");
        }
        {
            if(const auto* ffptr = dynamic_cast<ForceField<traitsT> const*>(ff.get()))
            {
                if(!ffptr->external().empty())
                {
                    MJOLNIR_LOG_WARN("GFW NPT langevin integrator currently does not "
                        "consider virial from external forcefield.");
                }
            }
            else
            {
                // MultipleBasin is too complicated and it is difficult to check
                // all the basins in all the units
                MJOLNIR_LOG_WARN("GFW NPT langevin integrator currently does not "
                    "consider virial from external forcefield.");
            }
        }

        // calculate parameters for each particles
        this->update(sys);

        // calculate force
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
        }
        sys.virial() = matrix33_type(0,0,0, 0,0,0, 0,0,0);
        ff->calc_force(sys);

        // calculate the current pressure using the force calculated here
        const auto h_cell = sys.boundary().width();
        this->P_ins_ = this->calc_pressure(sys, h_cell);

        sys.attribute("volume") = X(h_cell) * Y(h_cell) * Z(h_cell);

        return;
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

            for(std::size_t i=0; i<sys.size(); ++i)
            {
                sys.velocity(i) = this->S(sys.velocity(i),
                                          sys.force(i) * sys.rmass(i),
                                          v_over_h, coef_S);
                sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
            }
            sys.virial() = matrix33_type(0,0,0, 0,0,0, 0,0,0);
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

        sys.attribute("volume") = det_h;
        sys.boundary().set_boundary(cell_origin, cell_origin + h_cell);

        // --------------------------------------------------------------------
        // calc force

        // update neighbor list; reduce margin, reconstruct the list if needed
        const auto max_displacement = std::sqrt(max_displacement_sq_1) +
                                      std::sqrt(max_displacement_sq_2);
        ff->reduce_margin(2 * max_displacement, sys);

        ff->calc_force(sys);

        // --------------------------------------------------------------------
        // update particle velocities

        {
            const auto v_over_h = hadamard_product(-v_cell_, inv_h);
            const auto coef_S = make_coordinate<coordinate_type>(
                    std::exp(halfdt_ * X(v_over_h)),
                    std::exp(halfdt_ * Y(v_over_h)),
                    std::exp(halfdt_ * Z(v_over_h)));

            for(std::size_t i=0; i<sys.size(); ++i)
            {
                sys.velocity(i) = this->S(sys.velocity(i),
                                          sys.force(i) * sys.rmass(i),
                                          v_over_h, coef_S);
            }
        }

        // --------------------------------------------------------------------
        // calc current pressure
        // using just updated velocities and virial from `calc_force`

        this->P_ins_ = this->calc_pressure(sys, h_cell);

        P_diff = P_ins_ - P_ref_;

        // --------------------------------------------------------------------
        // update cell velocities

        v_cell_ += hadamard_product(halfdt_ * rm_cell_,
                   hadamard_product(det_h * inv_h, P_diff) - chi_kBT_ * inv_h);

        return t + dt_;
    }

    real_type                     delta_t()    const noexcept {return dt_;}
    std::vector<real_type> const& parameters() const noexcept {return gammas_;}

  private:

    real_type check_attribute_exists(const system_type& sys, const std::string& attr)
    {
        if(!sys.has_attribute(attr))
        {
            throw_exception<std::out_of_range>("mjolnir::GFWNPTLangevinIntegrator: "
                "It requires attribute `", attr, "`, but `", attr,
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

        MJOLNIR_LOG_INFO("1 / M_cell = ", rm_cell_);
        MJOLNIR_LOG_INFO("beta_cell  = ", noise_coeff_cell_);

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
        const auto det_h = math::X(h_cell) * math::Y(h_cell) * math::Z(h_cell);

        // instantaneous pressure
        const auto& vir = sys.virial();
        auto Pins = math::make_coordinate<coordinate_type>(vir(0,0), vir(1,1), vir(2,2));
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto m = sys.mass(i);
            const auto r = sys.position(i);
            const auto v = sys.velocity(i);
            const auto f = sys.force(i);

            // here the original form is a tensor product. but the diagonal term
            // of the tensor product (cij = ai * bj) becomes hadamard product
            // (cii = ai * bi).

            Pins += m * math::hadamard_product(v, v);
        }
        return Pins / det_h;
    }

    // When the off-diagonal terms vanish, the Su and the Sl functions become
    // the same with each other and are drastically simplified.
    // Here it implements diagonal version of Su and Sl.
    //
    // Since S always uses t = dt/2, here it assumes that t = dt/2.
    coordinate_type S(coordinate_type x0,      const coordinate_type b,
                      const coordinate_type A, const coordinate_type coef_S)
    {
        using math::X; using math::Y; using math::Z;

        X(x0) *= X(coef_S);
        Y(x0) *= Y(coef_S);
        Z(x0) *= Z(coef_S);

        // use tailor expansion of F1(0, adt/2), as in the paper
        const auto tA = halfdt_ * A;

        X(x0) += halfdt_ * X(b) * (real_type(1) + X(tA) / real_type(2) + X(tA) * X(tA) / real_type(6));
        Y(x0) += halfdt_ * Y(b) * (real_type(1) + Y(tA) / real_type(2) + Y(tA) * Y(tA) / real_type(6));
        Z(x0) += halfdt_ * Z(b) * (real_type(1) + Z(tA) / real_type(2) + Z(tA) * Z(tA) / real_type(6));

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
extern template class GFWNPTLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class GFWNPTLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class GFWNPTLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class GFWNPTLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif /* MJOLNIR_GFW_NPT_LANGEVIN_INTEGRATOR */
