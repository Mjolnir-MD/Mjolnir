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
    ~NaivePairCalculation() override = default;
    NaivePairCalculation(NaivePairCalculation const&)            = default;
    NaivePairCalculation& operator=(NaivePairCalculation const&) = default;
    NaivePairCalculation(NaivePairCalculation&&)                 = default;
    NaivePairCalculation& operator=(NaivePairCalculation&&)      = default;

    bool valid() const noexcept override {return true;}

    void initialize (const system_type& sys, const potential_type& pot) override
    {
        this->make(sys, pot);
        return;
    }

    void make  (const system_type& sys, const potential_type& pot) override;
    void update(const real_type, const system_type&,
                const potential_type&) override
    {
        return;
    }

    real_type cutoff() const noexcept override {return std::numeric_limits<real_type>::infinity();}
    real_type margin() const noexcept override {return std::numeric_limits<real_type>::infinity();}

    range_type partners(std::size_t i) const noexcept override {return neighbors_[i];}

  private:

    neighbor_list_type  neighbors_;
};

template<typename traitsT, typename potentialT>
void NaivePairCalculation<traitsT, potentialT>::make(
        const system_type&, const potential_type& pot)
{
    this->neighbors_.clear();

    const auto& participants = pot.participants();

    std::vector<neighbor_type> partners;
    for(std::size_t idx=0; idx<participants.size(); ++idx)
    {
        partners.clear();
        const auto   i = participants[idx];
        for(std::size_t jdx=idx+1; jdx<participants.size(); ++jdx)
        {
            const auto j = participants[jdx];
            if(pot.has_interaction(i, j)) // likely
            {
                partners.emplace_back(j, pot.prepare_params(i, j));
            }
        }
        this->neighbors_.add_list_for(i, partners.begin(), partners.end());
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
extern template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>,        DebyeHuckelPotential<double>>;
extern template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>,        DebyeHuckelPotential<float >>;
extern template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<double>>;
extern template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<float >>;

extern template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>,        ExcludedVolumePotential<double>>;
extern template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>,        ExcludedVolumePotential<float >>;
extern template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<double>>;
extern template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<float >>;

extern template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>,        LennardJonesPotential<double>>;
extern template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>,        LennardJonesPotential<float >>;
extern template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<double>>;
extern template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<float >>;

extern template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>,        UniformLennardJonesPotential<double>>;
extern template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>,        UniformLennardJonesPotential<float >>;
extern template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<double>>;
extern template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<float >>;
}
#endif // SEPARATE_BUILD

#endif /*MJOLNIR_CORE_NAIVE_PAIR_CALCULATION*/
