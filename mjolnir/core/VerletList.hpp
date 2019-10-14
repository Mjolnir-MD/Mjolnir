#ifndef MJOLNIR_CORE_VERLET_LIST_HPP
#define MJOLNIR_CORE_VERLET_LIST_HPP
#include <mjolnir/core/SpatialPartitionBase.hpp>
#include <algorithm>
#include <limits>

namespace mjolnir
{

template<typename traitsT, typename PotentialT, bool Newtons3rdLaw = true>
class VerletList final : public SpatialPartitionBase<traitsT, PotentialT>
{
  public:

    using traits_type        = traitsT;
    using potential_type     = PotentialT;
    using base_type          = SpatialPartitionBase<traits_type, potential_type>;

    using system_type        = typename base_type::system_type;
    using boundary_type      = typename base_type::boundary_type;
    using real_type          = typename base_type::real_type;
    using coordinate_type    = typename base_type::coordinate_type;
    using neighbor_list_type = typename base_type::neighbor_list_type;
    using neighbor_type      = typename base_type::neighbor_type;
    using range_type         = typename base_type::range_type;

  public:

    VerletList() : margin_(0.5), current_margin_(-1.0){}
    explicit VerletList(const real_type mgn): margin_(mgn), current_margin_(-1.0){}

    ~VerletList() override {}
    VerletList(VerletList const&) = default;
    VerletList(VerletList &&)     = default;
    VerletList& operator=(VerletList const&) = default;
    VerletList& operator=(VerletList &&)     = default;

    bool valid() const noexcept override
    {
        return current_margin_ >= 0.0;
    }

    void initialize(neighbor_list_type& neighbors,
            const system_type& sys, const potential_type& pot) override
    {
        this->set_cutoff(pot.max_cutoff_length());
        this->make(neighbors, sys, pot);
        return;
    }

    void make  (neighbor_list_type& neighbors,
                const system_type& sys, const potential_type& pot) override;

    void update(neighbor_list_type& neighbors, const real_type dmargin,
                const system_type& sys, const potential_type& pot) override
    {
        this->current_margin_ -= dmargin;
        if(this->current_margin_ < 0)
        {
            this->make(neighbors, sys, pot);
        }
        return ;
    }

    real_type cutoff() const noexcept override {return this->cutoff_;}
    real_type margin() const noexcept override {return this->margin_;}

  private:

    void set_cutoff(const real_type c) noexcept {this->cutoff_ = c;}
    void set_margin(const real_type m) noexcept {this->margin_ = m;}

  private:

    real_type      cutoff_;
    real_type      margin_;
    real_type      current_margin_;
};

template<typename traitsT, typename potentialT, bool N3L>
void VerletList<traitsT, potentialT, N3L>::make(neighbor_list_type& neighbors,
        const system_type& sys, const potential_type& pot)
{
    neighbors.clear();

    // `participants` is a list that contains indices of particles that are
    // related to the potential.
    const auto& participants = pot.participants();

    const real_type rc = cutoff_ * (1. + margin_);
    const real_type rc2 = rc * rc;

    std::vector<neighbor_type> partner;
    for(std::size_t idx=0; idx<participants.size(); ++idx)
    {
        partner.clear();
        const auto   i = participants[idx];
        const auto& ri = sys.position(i);

        // When N3L flag is activated, (j, i) is omitted if i < j. Since forces
        // on an interacting pair are always equal in magnitude and opposite in
        // direction, we normally don't need to calculate forces twice.
        //     In some cases, like a special forcefield, requires a complete
        // list of interacting pair.
        for(std::size_t jdx = (N3L ? idx+1 : 0); jdx<participants.size(); ++jdx)
        {
            const auto j = participants[jdx];
            if(!pot.has_interaction(i, j))
            {
                continue;
            }
            const auto& rj = sys.position(j);
            if(math::length_sq(sys.adjust_direction(rj - ri)) < rc2)
            {
                partner.emplace_back(j, pot.prepare_params(i, j));
            }
        }
        // because j is searched sequencially, sorting is not needed.
        neighbors.add_list_for(i, partner.begin(), partner.end());
    }
    this->current_margin_ = cutoff_ * margin_;
    return ;
}
} // mjolnir


#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/potential/global/DebyeHuckelPotential.hpp>
#include <mjolnir/potential/global/ExcludedVolumePotential.hpp>
#include <mjolnir/potential/global/LennardJonesPotential.hpp>
#include <mjolnir/potential/global/UniformLennardJonesPotential.hpp>

namespace mjolnir
{
extern template class VerletList<SimulatorTraits<double, UnlimitedBoundary>,        DebyeHuckelPotential<double>, true>;
extern template class VerletList<SimulatorTraits<float,  UnlimitedBoundary>,        DebyeHuckelPotential<float >, true>;
extern template class VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>, true>;
extern template class VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float >, true>;

extern template class VerletList<SimulatorTraits<double, UnlimitedBoundary>,        ExcludedVolumePotential<double>, true>;
extern template class VerletList<SimulatorTraits<float,  UnlimitedBoundary>,        ExcludedVolumePotential<float >, true>;
extern template class VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>, true>;
extern template class VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float >, true>;

extern template class VerletList<SimulatorTraits<double, UnlimitedBoundary>,        LennardJonesPotential<double>, true>;
extern template class VerletList<SimulatorTraits<float,  UnlimitedBoundary>,        LennardJonesPotential<float >, true>;
extern template class VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>, true>;
extern template class VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float >, true>;

extern template class VerletList<SimulatorTraits<double, UnlimitedBoundary>,        UniformLennardJonesPotential<double>, true>;
extern template class VerletList<SimulatorTraits<float,  UnlimitedBoundary>,        UniformLennardJonesPotential<float >, true>;
extern template class VerletList<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>, true>;
extern template class VerletList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float >, true>;
}
#endif // SEPARATE_BUILD

#endif/* MJOLNIR_CORE_VERLET_LIST */
