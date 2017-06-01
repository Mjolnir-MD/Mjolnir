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
    typedef std::vector<std::size_t> index_array;
    typedef std::vector<index_array> partners_type;

    struct information
    {
        information() : chain_idx(std::numeric_limits<std::size_t>::max()){}
        std::size_t chain_idx;
        index_array except_chains;
        index_array except_indices;
    };
    typedef std::vector<information> particle_info_type;

  public:

    NaivePairCalculation() = default;
    ~NaivePairCalculation() = default;

    NaivePairCalculation(NaivePairCalculation&&)            = default;
    NaivePairCalculation& operator=(NaivePairCalculation&&) = default;

    bool valid() const noexcept {return true;}

    std::size_t& chain_index   (std::size_t i);
    index_array& except_indices(std::size_t i);
    index_array& except_chains (std::size_t i);

    void initialize(const system_type& sys) const noexcept {return;}
    void make  (const system_type& sys);
    void update(const system_type& sys) noexcept {return;}
    void update(const system_type& sys, const real_type dt) noexcept {return;}

    index_array const& partners(std::size_t i) const noexcept {return partners_[i];}

  private:

    partners_type      partners_;
    particle_info_type informations_;
};

template<typename traitsT>
std::size_t& NaivePairCalculation<traitsT>::chain_index(std::size_t i)
{
    if(this->informations_.size() <= i)
        this->informations_.resize(i+1);
    return this->informations_.at(i).chain_idx;
}

template<typename traitsT>
typename NaivePairCalculation<traitsT>::index_array&
NaivePairCalculation<traitsT>::except_indices(std::size_t i)
{
    if(this->informations_.size() <= i)
        this->informations_.resize(i+1);
    return this->informations_.at(i).except_indices;
}

template<typename traitsT>
typename NaivePairCalculation<traitsT>::index_array&
NaivePairCalculation<traitsT>::except_chains(std::size_t i)
{
    if(this->informations_.size() <= i)
        this->informations_.resize(i+1);
    return this->informations_.at(i).except_chains;
}

template<typename traitsT>
void NaivePairCalculation<traitsT>::make(const system_type& sys)
{
    this->partners_.resize(sys.size());
    for(auto& partner : this->partners) partner.clear();

    if(informations_.size() < sys.size())
        informations_.resize(sys.size());

    for(std::size_t i=0, sz = sys.size()-1; i < sz; ++i)
    {
        const auto& info = informations_.at(i);
        const auto index_begin = info.except_indices.cbegin();
        const auto index_end   = info.except_indices.cend();
        const auto chain_begin = info.except_chains.cbegin();
        const auto chain_end   = info.except_chains.cend();

        for(std::size_t j=i+1; j<sys.size(); ++j)
        {
            const std::size_t j_chain = informations_.at(j).chain_idx;
            if(std::find(chain_begin, chain_end, j_chain) != chain_end) continue;
            if(std::find(index_begin, index_end, j)       != index_end) continue;
            this->partners_.at(i).push_back(j);
        }
    }
    return;
}

} // mjolnir
#endif /*MJOLNIR_CORE_NAIVE_PAIR_CALCULATION*/
