#ifndef MJOLNIR_INTERACTION_LOCAL_DIRECTIONAL_CONTACT_INTERACTION_HPP
#define MJOLNIR_INTERACTION_LOCAL_DIRECTIONAL_CONTACT_INTERACTION_HPP
#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/logger.hpp>
#include <tuple>
#include <cmath>


namespace mjolnir
{

template<typename traitsT,           typename angle1_potentialT,
         typename angle2_potentialT, typename contact_potentialT>
class DirectionalContactInteraction final : public LocalInteractionBase<traitsT>
{
  public:
    using traits_type            = traitsT;
    using angle1_potential_type  = angle1_potentialT;
    using angle2_potential_type  = angle2_potentialT;
    using contact_potential_type = contact_potentialT;
    using base_type              = LocalInteractionBase<traits_type>;
    using real_type              = typename base_type::real_type;
    using coordinate_type        = typename base_type::coordinate_type;
    using system_type            = typename base_type::system_type;
    using topology_type          = typename base_type::topology_type;
    using connection_kind_type   = typename base_type::connection_kind_type;

    using indices_type             = std::array<std::size_t, 4>;
    using indices_potentials_tuple = std::tuple<indices_type,
          angle1_potential_type, angle2_potential_type, contact_potential_type>;
    using container_type = std::vector<indices_potentials_tuple>;

  public:

    DirectionalContactInteraction(const connection_kind_type kind,
                                  const container_type&      pot,
                                  const real_type            margin = 0.5)
        : kind_(kind), potentials_(pot), margin_(margin)
    {}
    DirectionalContactInteraction(const connection_kind_type kind,
                                  container_type&&           pot,
                                  const real_type            margin = 0.5)
        : kind_(kind), potentials_(pot), margin_(margin)
    {}
    ~DirectionalContactInteraction() override{}

    void      calc_force (system_type&)       const noexcept override;
    real_type calc_energy(const system_type&) const noexcept override;

    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("angle1 = ",     angle1_potential_type::name(),
                         ", angle2 = ",   angle2_potential_type::name(),
                         ", contact = ",  contact_potential_type::name(),
                         ", number of contacts = ", potentials_.size());
        this->cutoff_ = this->max_cutoff_length();
        this->make_list(sys);
        for(auto& potential : this->potentials_)
        {
            std::get<1>(potential).initialize(sys);
            std::get<2>(potential).initialize(sys);
            std::get<3>(potential).initialize(sys);
        }
        return;
    }

    void update(const system_type& sys) override
    {
        for(auto& item: potentials_)
        {
            std::get<1>(item).update(sys);
            std::get<2>(item).update(sys);
            std::get<3>(item).update(sys);
        }
        this->cutoff_ = this->max_cutoff_length();
        this->make_list(sys);
        return;
    }

    void update_margin(const real_type dmargin, const system_type& sys) override
    {
        this->current_margin_ -= dmargin;
        if(this->current_margin_ < 0)
        {
            this->make_list(sys);
        }
        return;
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        this->current_margin_ = (cutoff_ + current_margin_) * scale - cutoff_;
        if(this->current_margin_ < 0)
        {
            this->make_list(sys);
        }
        return;
    }



    std::string name() const override
    {
        return "DirectionalContact:"_s + angle1_potential_type::name() + ","_s +
                                         angle2_potential_type::name() + ","_s +
                                         contact_potential_type::name();
    }

    // DirectionalContact is a Contact, so the topology is defined between
    // indices[1] and [2].
    void write_topology(topology_type& topol) const override
    {
        if(this->kind_.empty() || this->kind_ == "none") {return;}

        for(const auto& idxp : this->potentials_)
        {
            const auto& indices = std::get<0>(idxp);
            const auto  Pi = indices[1];
            const auto  Pj = indices[2];
            topol.add_connection(Pi, Pj, this->kind_);
        }
        return;
    }

    container_type const& potentials() const noexcept {return potentials_;}
    container_type&       potentials()       noexcept {return potentials_;}

    base_type* clone() const override
    {
        return new DirectionalContactInteraction(
                kind_, container_type(potentials_));
    }

  private:

    void make_list(const system_type& sys)
    {
        this->active_contacts_.clear();
        this->active_contacts_.reserve(potentials_.size());

        // absolute length of margin (this->margin_ is a relative length).
        const real_type abs_margin = this->cutoff_ * this->margin_;

        for(std::size_t i=0; i<this->potentials_.size(); ++i)
        {
            const auto& pot  = this->potentials_[i];
            const auto& pos0 = sys.position(std::get<0>(pot)[1]);
            const auto& pos1 = sys.position(std::get<0>(pot)[2]);
            const auto  dpos = sys.adjust_direction(pos1 - pos0);
            const auto  len2 = math::length_sq(dpos);

            const real_type rc = std::get<3>(pot).cutoff() + abs_margin;
            if(len2 < rc * rc)
            {
                this->active_contacts_.push_back(i);
            }
        }
        this->current_margin_ = abs_margin;
        return;
    }

    real_type max_cutoff_length() const noexcept
    {
        const auto max_cutoff_potential_itr = std::max_element(
            potentials_.begin(), potentials_.end(),
            [](const indices_potentials_tuple& lhs,
               const indices_potentials_tuple& rhs)
            {
                return std::get<3>(lhs).cutoff() < std::get<3>(rhs).cutoff();
            });
        return std::get<3>(*max_cutoff_potential_itr).cutoff();
    }

  private:

    connection_kind_type kind_;
    container_type potentials_;

    // neighbor list stuff
    real_type cutoff_;
    real_type margin_;
    real_type current_margin_;
    std::vector<std::size_t> active_contacts_;
};

template<typename traitsT,           typename angle1_potentialT,
         typename angle2_potentialT, typename contact_potentialT>
void DirectionalContactInteraction<
    traitsT, angle1_potentialT, angle2_potentialT, contact_potentialT
    >::calc_force(system_type& sys) const noexcept
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();

    for(const std::size_t active_contact : this->active_contacts_)
    {
        const auto& idxp = this->potentials_[active_contact];

        const auto& angle1_pot  = std::get<1>(idxp);
        const auto& angle2_pot  = std::get<2>(idxp);
        const auto& contact_pot = std::get<3>(idxp);

        const std::size_t      Ci  = std::get<0>(idxp)[0];
        const std::size_t      Pi  = std::get<0>(idxp)[1];
        const std::size_t      Pj  = std::get<0>(idxp)[2];
        const std::size_t      Cj  = std::get<0>(idxp)[3];
        const coordinate_type& rCi = sys.position(Ci);
        const coordinate_type& rPi = sys.position(Pi);
        const coordinate_type& rPj = sys.position(Pj);
        const coordinate_type& rCj = sys.position(Cj);

        // =========================================================
        // contact schema
        //
        //     theta1 theta2
        //       |      |
        //  Ci o v      v o Cj
        //     \-.      ,-/
        //   Pi o- - - - o Pj
        //        |Pij|
        //
        // U_dir = U_angle1(theta1) * U_angle2(theta2) * U_contact(|Pij|)

        const auto  Pij = sys.adjust_direction(rPj - rPi); // Pi -> Pj
        const auto lPij = math::length(Pij);
        if(lPij > contact_pot.cutoff())
        {
            continue;
        }

        constexpr auto abs_tol = math::abs_tolerance<real_type>();

        // ==========================================================
        // dU_angle1(theta1) / dr
        const auto PiCi         = sys.adjust_direction(rCi - rPi);
        const auto inv_len_PiCi = math::rlength(PiCi);
        const auto PiCi_reg     = PiCi * inv_len_PiCi;

        const auto inv_len_Pij  = real_type(1.0) / lPij;
        const auto Pij_reg      = Pij * inv_len_Pij;

        const auto PiCi_dot_Pij = math::dot_product(PiCi_reg, Pij_reg);
        const auto cos_theta1   = math::clamp<real_type>(PiCi_dot_Pij, -1, 1);
        const auto theta1       = std::acos(cos_theta1);
        const auto angle1_coef  = angle1_pot.derivative(theta1);

        const auto sin_theta1          = std::sin(theta1);
        const auto angle1_coef_inv_sin = angle1_coef / std::max(sin_theta1, abs_tol);

        const auto dU_angle1_drCi = (angle1_coef_inv_sin * inv_len_PiCi) *
                                    (cos_theta1 * PiCi_reg - Pij_reg);
        const auto dU_angle1_drPj = (angle1_coef_inv_sin * inv_len_Pij)  *
                                    (cos_theta1 * Pij_reg  - PiCi_reg);

        const auto dU_angle1_drPi = -(dU_angle1_drCi + dU_angle1_drPj);

        MJOLNIR_LOG_DEBUG("dU_angle1 / drCi = {", dU_angle1_drCi, "}");
        MJOLNIR_LOG_DEBUG("dU_angle1 / drPj = {", dU_angle1_drPj, "}");
        MJOLNIR_LOG_DEBUG("dU_angle1 / drPi = {", dU_angle1_drPi, "}");

        // dU_angle2(theta2) / dr
        const auto PjCj         = sys.adjust_direction(rCj - rPj);
        const auto inv_len_PjCj = math::rlength(PjCj);
        const auto PjCj_reg     = PjCj * inv_len_PjCj;

        const auto Pji_reg      = -Pij_reg;
        const auto PjCj_dot_Pji = math::dot_product(PjCj_reg, Pji_reg);
        const auto cos_theta2   = math::clamp<real_type>(PjCj_dot_Pji, -1, 1);
        const auto theta2       = std::acos(cos_theta2);
        const auto angle2_coef  = angle2_pot.derivative(theta2);

        const auto sin_theta2          = std::sin(theta2);
        const auto angle2_coef_inv_sin = angle2_coef / std::max(sin_theta2, abs_tol);

        const auto dU_angle2_drCj = (angle2_coef_inv_sin * inv_len_PjCj) *
                                    (cos_theta2 * PjCj_reg - Pji_reg);
        const auto dU_angle2_drPi = (angle2_coef_inv_sin * inv_len_Pij)  *
                                    (cos_theta2 * Pji_reg  - PjCj_reg);

        const auto dU_angle2_drPj = -(dU_angle2_drCj + dU_angle2_drPi);

        MJOLNIR_LOG_DEBUG("dU_angle2 / drCj = {", dU_angle2_drCj, "}");
        MJOLNIR_LOG_DEBUG("dU_angle2 / drPi = {", dU_angle2_drPi, "}");
        MJOLNIR_LOG_DEBUG("dU_angle2 / drPj = {", dU_angle2_drPj, "}");

        // dU_con(|Pij|) / dr
        const auto contact_coef = contact_pot.derivative(lPij);
        const auto dU_con_drPj  = contact_coef * Pij_reg;
        const auto dU_con_drPi  = -dU_con_drPj;

        MJOLNIR_LOG_DEBUG("dU_con / drPi = {", dU_con_drPi, "}");
        MJOLNIR_LOG_DEBUG("dU_con / drPj = {", dU_con_drPj, "}");

        const real_type U_angle1          =  angle1_pot.potential(theta1);
        const real_type U_angle2          =  angle2_pot.potential(theta2);
        const real_type U_con             = contact_pot.potential(lPij);
        const real_type U_angle1_U_con    = U_angle1 * U_con;
        const real_type U_angle2_U_con    = U_angle2 * U_con;
        const real_type U_angle1_U_angle2 = U_angle1 * U_angle2;

        MJOLNIR_LOG_DEBUG("U_angle1 = ", U_angle1);
        MJOLNIR_LOG_DEBUG("U_angle2 = ", U_angle2);
        MJOLNIR_LOG_DEBUG("U_con    = ", U_con);

        const auto dU_dir_drCi = dU_angle1_drCi * U_angle2_U_con;
        const auto dU_dir_drPi = dU_angle1_drPi * U_angle2_U_con +
                                 dU_angle2_drPi * U_angle1_U_con +
                                 U_angle1_U_angle2 * dU_con_drPi;
        const auto dU_dir_drPj = dU_angle1_drPj * U_angle2_U_con +
                                 dU_angle2_drPj * U_angle1_U_con +
                                 U_angle1_U_angle2 * dU_con_drPj;
        const auto dU_dir_drCj = dU_angle2_drCj * U_angle1_U_con;

        MJOLNIR_LOG_DEBUG("dU_dir / drCi = {", dU_dir_drCi, "}");
        MJOLNIR_LOG_DEBUG("dU_dir / drPi = {", dU_dir_drPi, "}");
        MJOLNIR_LOG_DEBUG("dU_dir / drPj = {", dU_dir_drPj, "}");
        MJOLNIR_LOG_DEBUG("dU_dir / drCj = {", dU_dir_drCj, "}");

        sys.force(Ci) -= dU_dir_drCi;
        sys.force(Pi) -= dU_dir_drPi;
        sys.force(Pj) -= dU_dir_drPj;
        sys.force(Cj) -= dU_dir_drCj;
    }
    return;
}

template<typename traitsT,           typename angle1_potentialT,
         typename angle2_potentialT, typename contact_potentialT>
typename DirectionalContactInteraction<traitsT, angle1_potentialT,
         angle2_potentialT, contact_potentialT>::real_type
DirectionalContactInteraction<
    traitsT, angle1_potentialT, angle2_potentialT, contact_potentialT
    >::calc_energy(const system_type& sys) const noexcept
{
    real_type E = 0.0;
    for(const std::size_t active_contact : active_contacts_)
    {
        const auto& idxp = this->potentials_[active_contact];

        const auto& angle1_pot  = std::get<1>(idxp);
        const auto& angle2_pot  = std::get<2>(idxp);
        const auto& contact_pot = std::get<3>(idxp);

        const std::size_t      Ci  = std::get<0>(idxp)[0];
        const std::size_t      Pi  = std::get<0>(idxp)[1];
        const std::size_t      Pj  = std::get<0>(idxp)[2];
        const std::size_t      Cj  = std::get<0>(idxp)[3];
        const coordinate_type& rCi = sys.position(Ci);
        const coordinate_type& rPi = sys.position(Pi);
        const coordinate_type& rPj = sys.position(Pj);
        const coordinate_type& rCj = sys.position(Cj);

        const auto  Pij = sys.adjust_direction(rPj - rPi); // Pi -> Pj
        const auto lPij = math::length(Pij);
        if(lPij > contact_pot.cutoff())
        {
            continue;
        }

        // calculate theta1
        const auto PiCi         = sys.adjust_direction(rCi - rPi);
        const auto inv_len_PiCi = math::rlength(PiCi);
        const auto PiCi_reg     = PiCi * inv_len_PiCi;

        const auto inv_len_Pij  = real_type(1.0) / lPij;
        const auto Pij_reg      = Pij * inv_len_Pij;

        const auto PiCi_dot_Pij = math::dot_product(PiCi_reg, Pij_reg);
        const auto cos_theta1   = math::clamp<real_type>(PiCi_dot_Pij, -1, 1);
        const auto theta1       = std::acos(cos_theta1);

        // calculate theta2
        const auto PjCj         = sys.adjust_direction(rCj - rPj);
        const auto inv_len_PjCj = math::rlength(PjCj);
        const auto PjCj_reg     = PjCj * inv_len_PjCj;

        const auto Pji_reg      = -Pij_reg;
        const auto PjCj_dot_Pji = math::dot_product(PjCj_reg, Pji_reg);
        const auto cos_theta2   = math::clamp<real_type>(PjCj_dot_Pji, -1, 1);
        const auto theta2       = std::acos(cos_theta2);

        E += angle1_pot.potential(theta1) * angle2_pot.potential(theta2) *
             contact_pot.potential(lPij);
    }
    return E;
}

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/potential/local/CosinePotential.hpp>
#include <mjolnir/potential/local/UniformPotential.hpp>
#include <mjolnir/potential/local/GoContactPotential.hpp>
#include <mjolnir/potential/local/GaussianPotential.hpp>

namespace mjolnir
{

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, CosinePotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , CosinePotential<float> , GaussianPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, CosinePotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , CosinePotential<float> , GaussianPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, GaussianPotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , GaussianPotential<float> , GaussianPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, GaussianPotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , GaussianPotential<float> , GaussianPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, UniformPotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , UniformPotential<float> , GaussianPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, UniformPotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , UniformPotential<float> , GaussianPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, GaussianPotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , GaussianPotential<float> , GaussianPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, GaussianPotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , GaussianPotential<float> , GaussianPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, UniformPotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , UniformPotential<float> , GaussianPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, UniformPotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , UniformPotential<float> , GaussianPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, CosinePotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , CosinePotential<float> , GaussianPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, CosinePotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , CosinePotential<float> , GaussianPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, UniformPotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , UniformPotential<float> , GaussianPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, UniformPotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , UniformPotential<float> , GaussianPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, CosinePotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , CosinePotential<float> , GaussianPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, CosinePotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , CosinePotential<float> , GaussianPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, GaussianPotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , GaussianPotential<float> , GaussianPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, GaussianPotential<double>, GaussianPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , GaussianPotential<float> , GaussianPotential<float> >;

// ---------------------------------------------------------

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, CosinePotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , CosinePotential<float> , GoContactPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, CosinePotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , CosinePotential<float> , GoContactPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, GaussianPotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , GaussianPotential<float> , GoContactPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, GaussianPotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , GaussianPotential<float> , GoContactPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, UniformPotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , UniformPotential<float> , GoContactPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, UniformPotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , UniformPotential<float> , GoContactPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, GaussianPotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , GaussianPotential<float> , GoContactPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, GaussianPotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , GaussianPotential<float> , GoContactPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, UniformPotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , UniformPotential<float> , GoContactPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, UniformPotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , UniformPotential<float> , GoContactPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, CosinePotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , CosinePotential<float> , GoContactPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, CosinePotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , CosinePotential<float> , GoContactPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, UniformPotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , UniformPotential<float> , GoContactPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, UniformPotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , UniformPotential<float> , GoContactPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, CosinePotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , CosinePotential<float> , GoContactPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, CosinePotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , CosinePotential<float> , GoContactPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, GaussianPotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , GaussianPotential<float> , GoContactPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, GaussianPotential<double>, GoContactPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , GaussianPotential<float> , GoContactPotential<float> >;

// --------------------------------------------

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, CosinePotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , CosinePotential<float> , UniformPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, CosinePotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , CosinePotential<float> , UniformPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, GaussianPotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , GaussianPotential<float> , UniformPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, GaussianPotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , GaussianPotential<float> , UniformPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, UniformPotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , UniformPotential<float> , UniformPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, UniformPotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , UniformPotential<float> , UniformPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, GaussianPotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , GaussianPotential<float> , UniformPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, GaussianPotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , GaussianPotential<float> , UniformPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, CosinePotential<double>, UniformPotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, CosinePotential<float> , UniformPotential<float> , UniformPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>, UniformPotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> , UniformPotential<float> , UniformPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, CosinePotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , CosinePotential<float> , UniformPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, CosinePotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , CosinePotential<float> , UniformPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GaussianPotential<double>, UniformPotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GaussianPotential<float> , UniformPotential<float> , UniformPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>, UniformPotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> , UniformPotential<float> , UniformPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, CosinePotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , CosinePotential<float> , UniformPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, CosinePotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , CosinePotential<float> , UniformPotential<float> >;

extern template class DirectionalContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, UniformPotential<double>, GaussianPotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, UniformPotential<float> , GaussianPotential<float> , UniformPotential<float> >;
extern template class DirectionalContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformPotential<double>, GaussianPotential<double>, UniformPotential<double>>;
extern template class DirectionalContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformPotential<float> , GaussianPotential<float> , UniformPotential<float> >;

} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_INTERACTION_LOCAL_DIRECTIONAL_CONTACT_INTERACTION_HPP */
