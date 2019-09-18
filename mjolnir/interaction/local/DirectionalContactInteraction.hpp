#ifndef MJOLNIR_INTERACTION_DIRECTIONAL_CONTACT_INTERACTION_HPP
#define MJOLNIR_INTERACTION_DIRECTIONAL_CONTACT_INTERACTION_HPP
#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/math/math.hpp>

namespace mjolnir
{

template<typename traitsT, typename angle1_potentialT, typename angle2_potentialT,
         typename contact_potentialT>
class DirectionalContactInteraction final : public LocalInteractionBase<traitsT>
{
 public:
  using traits_type          = traitsT;
  using angle1_pot_type      = angle1_potentialT;
  using angle2_pot_type      = angle2_potentialT;
  using contact_pot_type     = contact_potentialT;
  using base_type            = LocalInteractionBase<traits_type>;
  using real_type            = typename base_type::real_type;
  using coordinate_type      = typename base_type::coordinate_type;
  using system_type          = typename base_type::system_type;
  using topology_type        = typename base_type::topology_type;
  using connection_kind_type = typename base_type::connection_kind_type;

  using indices_type   = std::array<std::size_t, 4>;
  using indices_potentials_tuple = std::tuple<indices_type, angle1_pot_type,
                                              angle2_pot_type, contact_pot_type>;
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

  void      calc_force (system_type&) const noexcept override;
  real_type calc_energy(const system_type&) const noexcept override;

  void initialize(const system_type& sys) override
  {
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    MJOLNIR_LOG_INFO("angle1 potential = ", angle1_pot_type::name(),
                     ", angle2 potential = ", angle2_pot_type::name(),
                     ", contact potential = ", contact_pot_type::name(),
                     ", number of contacts = ", potentials_.size());

    this->cutoff_ = this->max_cutoff_length();
    this->make_list(sys);
    return;
  }

  void update(const system_type& sys) override
  {
    for(auto& item: potentials_)
    {
      std::get<contact_pot_type>(item).update(sys);
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

  std::string name() const override
  {
    return "DirectionalContact: {angle1 potentail: " + angle1_pot_type::name()
        + ", angle2 potential: " + angle2_pot_type::name()
        + ", contact potential: " + contact_pot_type::name() + "}";
  }

  void write_topology(topology_type&) const override;

 private:

  void make_list(const system_type& sys)
  {
    this->active_contacts_.clear();
    this->active_contacts_.reserve(potentials_.size());

    // absolute length of margin (this->margin_ is a relative length).
    const real_type abs_margin = this->cutoff_ * this->margin_;

    for(std::size_t i=0;i < this->potentials_.size(); ++i)
    {
      const indices_potentials_tuple& pot = this->potentials_[i];
      const coordinate_type pos0 = sys.position(std::get<indices_type>(pot).first[1]);
      const coordinate_type pos1 = sys.position(std::get<indices_type>(pot).first[2]);
      const coordinate_type dpos = sys.adjust_direction(pos1 - pos0);
      const real_type       len2 = math::length_sq(dpos);

      const real_type rc = std::get<contact_pot_type>(pot).cutoff() + abs_margin;
      if(len2 < rc * rc)
      {
        this->active_contacts_.push_back(i);
      }
      this->current_margin_ = this->cutoff_ * this->margin_;
      return;
    }
  }

  const real_type max_cutoff_length() const noexcept
  {
    const real_type max_cutoff = std::max_element(potentials_.begin(), potentials_.end(),
        [](const indices_potentials_tuple& lhs, const indices_potentials_tuple& rhs)
        {
          return lhs.second.cutoff() < rhs.second.cutoff();
        })->second.cutoff();
    return max_cutoff;
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

template<typename traitsT, typename angle1_potentialT,
         typename angle2_potentialT, typename contact_potentialT>
void DirectionalContactInteraction<traitsT, angle1_potentialT,
                                   angle2_potentialT, contact_potentialT>
::calc_force(system_type& sys) const noexcept
{
  for(const std::size_t active_contact : this -> active_contacts_)
  {
    const auto& idxp = this->potentials_[active_contact];

    const auto angle1_pot = std::get<angle1_pot_type>(idxp);
    const auto angle2_pot = std::get<angle2_pot_type>(idxp);
    const auto contact_pot = std::get<contact_pot_type>(idxp);

    const std::size_t Ci  = std::get<indices_type>(idxp)[0];
    const coordinate_type& rCi = sys.position(Ci);
    const std::size_t Pi = std::get<indices_type>(idxp)[1];
    const coordinate_type& rPi = sys.position(Pi);
    const std::size_t Cj  = std::get<indices_type>(idxp)[2];
    const coordinate_type& rCj = sys.position(Cj);
    const std::size_t Pj = std::get<indices_type>(idxp)[3];
    const coordinate_type& rPj = sys.position(Pj);

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

    const coordinate_type Pij = sys.adjust_direction(rPi - rPj); // Pi -> Pj
    const real_type lPij = math::length(Pij);
    if(lPij > contact_pot.cutoff_())
    {
      continue;
    }

    // ==========================================================
    // dU_angle1(theta1) / dr
    const coordinate_type PiCi = sys.adjust_direction(rCi - rPi);
    const real_type inv_len_PiCi = math::rlength(PiCi);
    const coordinate_type PiCi_reg = PiCi * inv_len_PiCi;

    const real_type inv_len_Pij = real_type(1.0) / lPij;
    const coordinate_type Pij_reg = Pij * inv_len_Pij;

    const real_type PiCi_dot_Pij = math::dot_product(PiCi_reg, Pij_reg);
    const real_type cos_theta1 = math::clamp(PiCi_dot_Pij,  real_type(-1.0), real_type(1.0));
    const real_type theta1 = std::acos(cos_theta1);
    const real_type angle1_coef = angle1_pot.derivative(theta1);

    const real_type sin_theta1 = std::sin(theta1);
    const real_type angle1_coef_inv_sin = (sin_theta1 > math::abs_tolerance<real_type>()) ?
        angle1_coef / sin_theta1 : angle1_coef / math::abs_tolerance<real_type>();

    const coordinate_type dU_angle1_drCi =
        (angle1_coef_inv_sin * inv_len_PiCi) * (cos_theta1 * PiCi_reg - Pij_reg);
    const coordinate_type dU_angle1_drPj =
        (angle1_coef_inv_sin * inv_len_Pij) * (cos_theta1 * Pij_reg - PiCj_reg);

    const coordinate_type dU_angle1_drPi = -(dU_angle1_drCi + dU_angle1_drPj);

    // dU_angle2(theta2) / dr
    const coordinate_type PjCj = sys.adjust_direction(rPj - rCj);
    const real_type inv_len_PjCj = math::rlength(PjCj);
    const coordinate_type PjCj_reg = PjCj * inv_len_PiCi;

    const Pji_reg = - Pij_reg;
    const real_type PjCj_dot_Pji = math::dot_product(PjCj_reg, Pji_reg);
    const real_type cos_theta2 = math::clamp(PjCj_dot_Pji, real_type(-1.0), real_type(1.0));
    const real_type theta2 = std::acos(cos_theta2);
    const real_type angle2_coef = angle2_pot.derivative(theta2);

    const real_type sin_theta2 = std::sin(theta2);
    const real_type angle2_coef_inv_sin = (sin_theta2 > math::abs_tolerance<real_type>()) ?
        angle2_coef / sin_theta2 : angle2_coef / math::abs_tolerance<real_type>();

    const coordinate_type dU_angle2_drCj =
        (angle2_coef_inv_sin * inv_len_PjCj) * (cos_theta2 * PjCj_reg - Pji_reg);
    const coordinate_type dU_angle2_drPi =
        (angle2_coef_inv_sin * inv_len_Pij) * (cos_theta2 * Pji_reg - PjCj_reg);

    const coordinate_type dU_angle2_drPj = -(dU_angle2_drCj + dU_angle2_drPj);

    // dU_con(|Pij|) / dr
    const real_type contact_coef = contact_pot.derivative(lPij);
    const coordinate_type dU_con_drPj = contact_coef * Pij_reg;
    const coordinate_type dU_con_drPi = - dU_con_drPj;

    const real_type U_angle1 = angle1_pot.calc_energy(theta1);
    const real_type U_angle2 = angle2_pot.calc_energy(theta2);
    const real_type U_con    = contact_pot.calc_energy(lPij);
    const real_type U_angle1_U_con = U_angle1 * U_con;
    const real_type U_angle2_U_con = U_angle2 * U_con;
    const real_type U_angle1_U_angle2 = U_angle1 * U_angle2;

    const coordinate_type force_Ci = dU_angle1_drCi * U_angle2_U_con;
    const coordinate_type force_Pi = dU_angle1_drPi * U_angle2_U_con
        + dU_angle2_drPi * U_angle1_U_con + U_angle1_U_angle2 * dU_con_drPi;
    const coordinate_type force_Pj = dU_angle1_drPj * U_angle2_U_con
        + dU_angle2_drPj * U_angle1_U_con + U_angle1_U_angle2 * dU_con_drPj;
    const coordinate_type force_Cj = dU_angle2_drCj * U_angle1_U_con;

    sys.force(Ci) -= force_Ci;
    sys.force(Pi) -= force_Pi;
    sys.force(Pj) -= force_Pj;
    sys.force(Cj) -= force_Cj;
  }
  return;
}

template<typename traitsT, typename angle1_potentialT,
         typename angle2_potentialT, typename contact_potentialT>
typename DirectionalContactInteraction<traitsT, angle1_potentialT,
                                       angle2_potentialT, contact_potentialT>::real_type
DirectionalContactInteraction<traitsT, angle1_potentialT, angle2_potentialT,
    contact_potentialT>::calc_energy(const system_type& sys) const noexcept
{
  real_type E = 0.0;
  return E;
}

template<typename traitsT, typename angle1_potentialT,
         typename angle2_potentialT, typename contact_potentialT>
void DirectionalContactInteraction<traitsT, angle1_potentialT,
                                   angle2_potentialT, contact_potentialT>
::write_topology(topology_type& topol) const
{
  return;
}

} // mjolnir

#endif /* MJOLNIR_DIRECTIONAL_CONTACT_INTERACTION */
