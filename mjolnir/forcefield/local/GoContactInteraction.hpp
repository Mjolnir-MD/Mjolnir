#ifndef MJOLNIR_INTERACTION_GO_CONTACT_INTERACTION_HPP
#define MJOLNIR_INTERACTION_GO_CONTACT_INTERACTION_HPP
#include <mjolnir/forcefield/local/ContactInteraction.hpp>
#include <mjolnir/forcefield/local/GoContactPotential.hpp>

namespace mjolnir
{

// It is a specialization of ContactInteraction for GoContactPotential.
// In the case of GoContactPotential, we can omit `sqrt` call that is
// normally used to calculate distance because we only needs the squared distance.
template<typename realT, template<typename, typename> class boundaryT>
class ContactInteraction<
    SimulatorTraits<realT, boundaryT>,
    GoContactPotential<realT>
    > final : public LocalInteractionBase<SimulatorTraits<realT, boundaryT>>
{
  public:
    using traits_type          = SimulatorTraits<realT, boundaryT>;
    using potential_type       = GoContactPotential<realT>;
    using base_type            = LocalInteractionBase<traits_type>;
    using real_type            = typename base_type::real_type;
    using coordinate_type      = typename base_type::coordinate_type;
    using system_type          = typename base_type::system_type;
    using topology_type        = typename base_type::topology_type;
    using connection_kind_type = typename base_type::connection_kind_type;

    using indices_type         = std::array<std::size_t, 2>;
    using potential_index_pair = std::pair<indices_type, potential_type>;
    using container_type       = std::vector<potential_index_pair>;
    using iterator             = typename container_type::iterator;
    using const_iterator       = typename container_type::const_iterator;

  public:

    ContactInteraction(const connection_kind_type kind,
                       const container_type&      pot,
                       const real_type            margin = 0.5)
        : kind_(kind), potentials_(pot), margin_(margin)
    {}
    ContactInteraction(const connection_kind_type kind,
                       container_type&&           pot,
                       const real_type            margin = 0.5)
        : kind_(kind), potentials_(std::move(pot)), margin_(margin)
    {}
    ~ContactInteraction() override {}

    real_type calc_energy(const system_type& sys) const noexcept override
    {
        real_type E = 0.0;
        for(const std::size_t active_contact : active_contacts_)
        {
            const auto& idxp = this->potentials_[active_contact];
            E += idxp.second.potential(math::length(sys.adjust_direction(
                    sys.position(idxp.first[0]), sys.position(idxp.first[1]))));
        }
        return E;
    }
    void calc_force(system_type& sys) const noexcept override
    {
        this->template calc_force_energy_virial_impl<false, false>(sys);
        return;
    }
    void calc_force_and_virial(system_type& sys) const noexcept override
    {
        this->template calc_force_energy_virial_impl<false, true>(sys);
        return;
    }
    real_type calc_force_and_energy(system_type& sys) const noexcept override
    {
        return this->template calc_force_energy_virial_impl<true, false>(sys);
    }
    real_type calc_force_virial_energy(system_type& sys) const noexcept override
    {
        return this->template calc_force_energy_virial_impl<true, true>(sys);
    }

    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential = ", potential_type::name(),
                         ", number of bonds = ", potentials_.size());

        this->cutoff_ = std::max_element(potentials_.begin(), potentials_.end(),
            [](const potential_index_pair& lhs, const potential_index_pair& rhs)
            {
                return lhs.second.cutoff() < rhs.second.cutoff();
            })->second.cutoff();
        this->make_list(sys);
        for(auto& potential : potentials_)
        {
            potential.second.initialize(sys);
        }
        return;
    }

    void update(const system_type& sys) override
    {
        for(auto& item : potentials_)
        {
            item.second.update(sys);
        }

        this->cutoff_ = std::max_element(potentials_.begin(), potentials_.end(),
            [](const potential_index_pair& lhs, const potential_index_pair& rhs)
            {
                return lhs.second.cutoff() < rhs.second.cutoff();
            })->second.cutoff();
        this->make_list(sys);
        return;
    }

    void reduce_margin(const real_type dmargin, const system_type& sys) override
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
    {return "Contact:"_s + potential_type::name();}

    void write_topology(topology_type& topol) const override
    {
        if(this->kind_.empty() || this->kind_ == "none") {return;}

        for(const auto& idxp : this->potentials_)
        {
            const auto i = idxp.first[0];
            const auto j = idxp.first[1];
            topol.add_connection(i, j, this->kind_);
        }
        return;
    }

    container_type const& potentials() const noexcept {return potentials_;}
    container_type&       potentials()       noexcept {return potentials_;}

    base_type* clone() const override
    {
        return new ContactInteraction(kind_, container_type(potentials_));
    }

  private:

    void make_list(const system_type& sys)
    {
        this->active_contacts_.clear();
        this->active_contacts_.reserve(potentials_.size());

        // absolute length of margin (this->margin_ is a relative length).
        const real_type abs_margin = this->cutoff_ * this->margin_;

        for(std::size_t i=0; i < this->potentials_.size(); ++i)
        {
            const auto& pot = this->potentials_[i];
            const auto pos0 = sys.position(pot.first[0]);
            const auto pos1 = sys.position(pot.first[1]);
            const auto dpos = sys.adjust_direction(pos0, pos1);
            const auto len2 = math::length_sq(dpos);

            const auto rc = pot.second.cutoff() + abs_margin;
            if(len2 < rc * rc)
            {
                this->active_contacts_.push_back(i);
            }
        }
        this->current_margin_ = this->cutoff_ * this->margin_;
        return;
    }

    template<bool NeedEnergy, bool NeedVirial>
    real_type calc_force_energy_virial_impl(system_type& sys) const noexcept
    {
        real_type energy = 0;
        for(const std::size_t active_contact : active_contacts_)
        {
            const auto& idxp = this->potentials_[active_contact];

            const std::size_t idx0 = idxp.first[0];
            const std::size_t idx1 = idxp.first[1];
            const auto&       pot  = idxp.second;

            const auto dpos =
                sys.adjust_direction(sys.position(idx0), sys.position(idx1));

            const real_type len2  = math::length_sq(dpos);
            if(pot.cutoff() * pot.cutoff() <= len2)
            {
                continue;
            }

            const real_type r2     = real_type(1) / len2;
            const real_type v0r_2  = pot.v0() * pot.v0() * r2;
            const real_type v0r_6  = v0r_2 * v0r_2 * v0r_2;
            const real_type v0r_10 = v0r_6 * v0r_2 * v0r_2;
            const real_type v0r_12 = v0r_10 * v0r_2;

            if(NeedEnergy)
            {
                energy += pot.k() * (5 * v0r_12 - 6 * v0r_10);
            }

            const auto coef = -60 * pot.k() * r2 * (v0r_10 - v0r_12);
            const auto f    = coef * dpos;
            sys.force(idx0) -= f;
            sys.force(idx1) += f;

            if(NeedVirial)
            {
                sys.virial() += math::tensor_product(dpos, f);
            }
        }
        return energy;
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

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

namespace mjolnir
{

// go-contact
extern template class ContactInteraction<SimulatorTraits<double, UnlimitedBoundary       >, GoContactPotential<double>>;
extern template class ContactInteraction<SimulatorTraits<float,  UnlimitedBoundary       >, GoContactPotential<float> >;
extern template class ContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GoContactPotential<double>>;
extern template class ContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GoContactPotential<float> >;

} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_BOND_LENGTH_INTERACTION */
