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
