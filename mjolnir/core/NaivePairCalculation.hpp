#ifndef MJOLNIR_CORE_NAIVE_PAIR_CALCULATION_HPP
#define MJOLNIR_CORE_NAIVE_PAIR_CALCULATION_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/NeighborList.hpp>
#include <mjolnir/util/empty.hpp>

namespace mjolnir
{

template<typename traitsT, typename parameterT>
class NaivePairCalculation
{
  public:
    using traits_type         = traitsT;
    using system_type         = System<traits_type>;
    using boundary_type       = typename traits_type::boundary_type;
    using real_type           = typename traits_type::real_type;
    using coordinate_type     = typename traits_type::coordinate_type;
    using parameter_type      = parameterT;
    using neighbor_list_type  = NeighborList<parameter_type>;
    using neighbor_type       = typename neighbor_list_type::neighbor_type;
    using range_type          = typename neighbor_list_type::range_type;

  public:

    NaivePairCalculation() = default;
    ~NaivePairCalculation() = default;
    NaivePairCalculation(NaivePairCalculation const&)            = default;
    NaivePairCalculation& operator=(NaivePairCalculation const&) = default;
    NaivePairCalculation(NaivePairCalculation&&)                 = default;
    NaivePairCalculation& operator=(NaivePairCalculation&&)      = default;

    bool valid() const noexcept {return true;}

    template<typename PotentialT>
    void initialize (const system_type& sys, const PotentialT& pot);

    template<typename PotentialT>
    void make  (const system_type& sys, const PotentialT& pot);

    template<typename PotentialT>
    void update(const real_type, const system_type&, const PotentialT&) noexcept
    {return;}

    range_type partners(std::size_t i) const noexcept {return neighbors_[i];}

  private:

    neighbor_list_type  neighbors_;
};

template<typename traitsT, typename parameterT>
template<typename PotentialT>
void NaivePairCalculation<traitsT, parameterT>::initialize(
        const system_type& sys, const PotentialT& pot)
{
    this->make(sys, pot);
    return;
}

template<typename traitsT, typename parameterT>
template<typename PotentialT>
void NaivePairCalculation<traitsT, parameterT>::make(
        const system_type&, const PotentialT& pot)
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

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>, empty_t>;
extern template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>, empty_t>;
extern template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, empty_t>;
extern template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, empty_t>;

extern template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>, double>;
extern template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>, float >;
extern template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, double>;
extern template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, float >;

extern template class NaivePairCalculation<SimulatorTraits<double, UnlimitedBoundary>, std::pair<double, double>>;
extern template class NaivePairCalculation<SimulatorTraits<float,  UnlimitedBoundary>, std::pair<float , float >>;
extern template class NaivePairCalculation<SimulatorTraits<double, CuboidalPeriodicBoundary>, std::pair<double, double>>;
extern template class NaivePairCalculation<SimulatorTraits<float,  CuboidalPeriodicBoundary>, std::pair<float , float >>;
#endif

} // mjolnir
#endif /*MJOLNIR_CORE_NAIVE_PAIR_CALCULATION*/
