#ifndef MJOLNIR_CORE_NAIVE_PAIR_CALCULATION_HPP
#define MJOLNIR_CORE_NAIVE_PAIR_CALCULATION_HPP
#include <mjolnir/core/SpatialPartitionBase.hpp>
#include <mjolnir/util/empty.hpp>

namespace mjolnir
{

template<typename traitsT, typename PotentialT>
class NaivePairCalculation final : public SpatialPartitionBase<traitsT, PotentialT>
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

    NaivePairCalculation() = default;
    ~NaivePairCalculation() override {}
    NaivePairCalculation(NaivePairCalculation const&)            = default;
    NaivePairCalculation& operator=(NaivePairCalculation const&) = default;
    NaivePairCalculation(NaivePairCalculation&&)                 = default;
    NaivePairCalculation& operator=(NaivePairCalculation&&)      = default;

    bool valid() const noexcept override {return true;}

    void initialize(neighbor_list_type& neighbor,
            const system_type& sys, const potential_type& pot) override
    {
        this->make(neighbor, sys, pot);
        return;
    }

    void make  (neighbor_list_type& neighbor,
                const system_type& sys, const potential_type& pot) override;
    void update(neighbor_list_type&, const real_type,
                const system_type&, const potential_type&) override
    {
        return;
    }
    void scale_margin(neighbor_list_type&, const real_type,
                      const system_type&, const potential_type&) override
    {
        return;
    }

    real_type cutoff() const noexcept override {return std::numeric_limits<real_type>::infinity();}
    real_type margin() const noexcept override {return std::numeric_limits<real_type>::infinity();}

    base_type* clone() const override {return new NaivePairCalculation();}
};

template<typename traitsT, typename potentialT>
void NaivePairCalculation<traitsT, potentialT>::make(neighbor_list_type& neighbors,
        const system_type&, const potential_type& pot)
{
    neighbors.clear();

    const auto leading_participants = pot.leading_participants();

    std::vector<neighbor_type> partners;
    for(std::size_t idx=0; idx<leading_participants.size(); ++idx)
    {
        partners.clear();
        const auto i = leading_participants[idx];
        for(const auto j : pot.possible_partners_of(idx, i))
        {
            if(pot.has_interaction(i, j)) // likely
            {
                partners.emplace_back(j, pot.prepare_params(i, j));
            }
        }
        neighbors.add_list_for(i, partners.begin(), partners.end());
    }
    return;
}

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/potential/global/DebyeHuckelPotential.hpp>
#include <mjolnir/potential/global/ExcludedVolumePotential.hpp>
#include <mjolnir/potential/global/LennardJonesPotential.hpp>
#include <mjolnir/potential/global/UniformLennardJonesPotential.hpp>

namespace mjolnir
{
extern template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>       , DebyeHuckelPotential<SimulatorTraits<double, UnlimitedBoundary>       >>;
extern template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>       , DebyeHuckelPotential<SimulatorTraits<float,  UnlimitedBoundary>       >>;
extern template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>       , ExcludedVolumePotential<SimulatorTraits<double, UnlimitedBoundary>       >>;
extern template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>       , ExcludedVolumePotential<SimulatorTraits<float,  UnlimitedBoundary>       >>;
extern template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>       , LennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>       >>;
extern template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>       , LennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>       >>;
extern template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>       , UniformLennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>       >>;
extern template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>       , UniformLennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>       >>;
extern template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
}
#endif // SEPARATE_BUILD

#endif /*MJOLNIR_CORE_NAIVE_PAIR_CALCULATION*/
