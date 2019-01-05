#ifndef MJOLNIR_CORE_NAIVE_PAIR_CALCULATION
#define MJOLNIR_CORE_NAIVE_PAIR_CALCULATION
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/NeighborList.hpp>
#include <mjolnir/core/ExclusionList.hpp>

namespace mjolnir
{

template<typename traitsT, typename parameterT>
class NaivePairCalculation
{
  public:

    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::boundary_type   boundary_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

    typedef ExclusionList exclusion_list_type;

    typedef parameterT parameter_type;
    typedef NeighborList<parameter_type> neighbor_list_type;
    typedef typename neighbor_list_type::neighbor_type neighbor_type;
    typedef typename neighbor_list_type::range_type    range_type;

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
    void reconstruct(const system_type& sys, const PotentialT& pot)
    {
        this->initialize(sys, pot); // do the same thing as `initialize`
        return;
    }

    template<typename PotentialT>
    void make  (const system_type& sys, const PotentialT& pot);

    template<typename PotentialT>
    void update(const real_type, const system_type&, const PotentialT&) noexcept
    {return;}

    range_type partners(std::size_t i) const noexcept {return neighbors_[i];}

  private:

    exclusion_list_type exclusion_;
    neighbor_list_type  neighbors_;
};

template<typename traitsT, typename parameterT>
template<typename PotentialT>
void NaivePairCalculation<traitsT, parameterT>::initialize(
        const system_type& sys, const PotentialT& pot)
{
    this->exclusion_.make(sys, pot);
    this->make(sys, pot);
    return;
}

template<typename traitsT, typename parameterT>
template<typename PotentialT>
void NaivePairCalculation<traitsT, parameterT>::make(
        const system_type& sys, const PotentialT& pot)
{
    this->neighbors_.clear();
    for(std::size_t i=0, sz = sys.size()-1; i < sz; ++i)
    {
        std::vector<neighbor_type> partners;
        for(std::size_t j=i+1; j<sys.size(); ++j)
        {
            if(this->exclusion_.is_excluded(i, j))
            {
                continue;
            }
            partners.emplace_back(j, pot.prepair_params(i, j));
        }
        this->neighbors_.add_list_for(i, partners.begin(), partners.end());
    }
    return;
}

} // mjolnir
#endif /*MJOLNIR_CORE_NAIVE_PAIR_CALCULATION*/
