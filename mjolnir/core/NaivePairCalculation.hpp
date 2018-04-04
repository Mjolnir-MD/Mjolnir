#ifndef MJOLNIR_CORE_NAIVE_PAIR_CALCULATION
#define MJOLNIR_CORE_NAIVE_PAIR_CALCULATION
#include "System.hpp"

namespace mjolnir
{

template<typename traitsT>
class NaivePairCalculation
{
  public:

    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::boundary_type   boundary_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

    typedef ExclusionList exclusion_list_type;
    typedef NeighborList  neighbor_list_type;
    typedef neighbor_list_type::range_type range_type;

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

    void make  (const system_type& sys);
    void update(const system_type& sys) noexcept {return;}

    // it does not have cutoff stuff
    void set_cutoff(const real_type) noexcept {return;}
    void set_mergin(const real_type) noexcept {return;}

    range_type partners(std::size_t i) const noexcept {return neighbors_[i];}

  private:

    exclusion_list_type exclusion_;
    neighbor_list_type  neighbors_;
};

template<typename traitsT>
template<typename PotentialT>
void NaivePairCalculation<traitsT>::initialize(
        const system_type& sys, const PotentialT& pot)
{
    this->exclusion_.make(sys, pot);
    return;
}

template<typename traitsT>
void NaivePairCalculation<traitsT>::make(const system_type& sys)
{
    this->partners_.resize(sys.size());
    for(auto& partner : this->partners_)
    {
        partner.clear();
    }

    for(std::size_t i=0, sz = sys.size()-1; i < sz; ++i)
    {
        for(std::size_t j=i+1; j<sys.size(); ++j)
        {
            if(this->exclusion_.is_excluded(i, j))
            {
                continue;
            }
            this->partners_[i].push_back(j);
        }
    }
    return;
}

} // mjolnir
#endif /*MJOLNIR_CORE_NAIVE_PAIR_CALCULATION*/
