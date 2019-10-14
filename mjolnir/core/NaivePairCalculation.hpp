#ifndef MJOLNIR_CORE_NAIVE_PAIR_CALCULATION_HPP
#define MJOLNIR_CORE_NAIVE_PAIR_CALCULATION_HPP
#include <mjolnir/core/SpatialPartitionBase.hpp>
#include <mjolnir/util/empty.hpp>

namespace mjolnir
{

// It is the most trivial and slowest SpatialPartition.
// It has *all* the possible interacting pairs, regardless of the distances.
// For debugging, it is useful to compare the results with this implementation.
template<typename traitsT, typename PotentialT, bool Newtons3rdLaw = true>
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

    real_type cutoff() const noexcept override {return std::numeric_limits<real_type>::infinity();}
    real_type margin() const noexcept override {return std::numeric_limits<real_type>::infinity();}
};

template<typename traitsT, typename potentialT, bool N3L>
void NaivePairCalculation<traitsT, potentialT, N3L>::make(
        neighbor_list_type& neighbors,
        const system_type&, const potential_type& pot)
{
    neighbors.clear();
    const auto& participants = pot.participants();

    std::vector<neighbor_type> partners;
    for(std::size_t idx=0; idx<participants.size(); ++idx)
    {
        partners.clear();
        const auto i = participants[idx];

        // When N3L flag is activated, (j, i) is omitted if i < j. Since forces
        // on an interacting pair are always equal in magnitude and opposite in
        // direction, we normally don't need to calculate forces twice.
        //     In some cases, like a special forcefield, requires a complete
        // list of interacting pair.
        for(std::size_t jdx=(N3L ? idx+1 : 0); jdx<participants.size(); ++jdx)
        {
            const auto j = participants[jdx];
            if(pot.has_interaction(i, j))
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
extern template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>,        DebyeHuckelPotential<double>, true>;
extern template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>,        DebyeHuckelPotential<float >, true>;
extern template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>, true>;
extern template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float >, true>;

extern template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>,        ExcludedVolumePotential<double>, true>;
extern template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>,        ExcludedVolumePotential<float >, true>;
extern template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>, true>;
extern template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float >, true>;

extern template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>,        LennardJonesPotential<double>, true>;
extern template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>,        LennardJonesPotential<float >, true>;
extern template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>, true>;
extern template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float >, true>;

extern template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>,        UniformLennardJonesPotential<double>, true>;
extern template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>,        UniformLennardJonesPotential<float >, true>;
extern template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>, true>;
extern template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float >, true>;
}
#endif // SEPARATE_BUILD

#endif /*MJOLNIR_CORE_NAIVE_PAIR_CALCULATION*/
