#ifndef MJOLNIR_CORE_VERLET_LIST
#define MJOLNIR_CORE_VERLET_LIST
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/NeighborList.hpp>
#include <algorithm>
#include <limits>

namespace mjolnir
{

template<typename traitsT>
class VerletList
{
  public:

    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::boundary_type   boundary_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

    typedef NeighborList           neighbor_list_type;
    typedef neighbor_list_type::range_type range_type;

    struct information
    {
        information() : chain_idx(std::numeric_limits<std::size_t>::max()){}
        std::size_t chain_idx;
        std::vector<std::size_t> except_chains;
        std::vector<std::size_t> except_indices;
    };
    typedef std::vector<information> particle_info_type;

  public:

    VerletList() = default;
    ~VerletList() = default;
    VerletList(VerletList const&) = default;
    VerletList(VerletList &&)     = default;
    VerletList& operator=(VerletList const&) = default;
    VerletList& operator=(VerletList &&)     = default;

    VerletList(const real_type cutoff, const real_type mergin)
        : dt_(0.), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.)
    {}
    VerletList(const real_type cutoff, const real_type mergin,
               const real_type dt)
        : dt_(dt), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.)
    {}

    bool valid() const noexcept
    {
        return current_mergin_ >= 0. || dt_ == 0.;
    }

    void initialize(const system_type&) const noexcept {return;}
    void make  (const system_type& sys);
    void update(const system_type& sys);
    void update(const system_type& sys, const real_type dt);

    std::size_t&              chain_index   (std::size_t i);
    std::vector<std::size_t>& except_indices(std::size_t i);
    std::vector<std::size_t>& except_chains (std::size_t i);

    real_type cutoff() const noexcept {return this->cutoff_;}
    real_type mergin() const noexcept {return this->mergin_;}

    void set_cutoff(const real_type c) noexcept {return this->cutoff_ = c;}
    void set_mergin(const real_type m) noexcept {return this->mergin_ = m;}

    range_type partners(std::size_t i) const noexcept {return neighbors_[i];}

  private:

    real_type      dt_;
    real_type      cutoff_;
    real_type      mergin_;
    real_type      current_mergin_;
    static Logger& logger_;

    neighbor_list_type neighbors_;
    particle_info_type informations_;
};

template<typename traitsT>
Logger& VerletList<traitsT>::logger_ = LoggerManager<char>::get_logger("VerletList");

template<typename traitsT>
std::size_t& VerletList<traitsT>::chain_index(std::size_t i)
{
    if(this->informations_.size() <= i)
        this->informations_.resize(i+1);
    return this->informations_.at(i).chain_idx;
}

template<typename traitsT>
std::vector<std::size_t>&
VerletList<traitsT>::except_indices(std::size_t i)
{
    if(this->informations_.size() <= i)
        this->informations_.resize(i+1);
    return this->informations_.at(i).except_indices;
}

template<typename traitsT>
std::vector<std::size_t>&
VerletList<traitsT>::except_chains(std::size_t i)
{
    if(this->informations_.size() <= i)
        this->informations_.resize(i+1);
    return this->informations_.at(i).except_chains;
}

template<typename traitsT>
void VerletList<traitsT>::make(const system_type& sys)
{
    this->neighbors_.clear();

    if(informations_.size() < sys.size()){informations_.resize(sys.size());}

    std::vector<std::size_t> tmp;

    const real_type rc = cutoff_ * (1. + mergin_);
    const real_type rc2 = rc * rc;
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto& ri = sys[i].position;
        const auto& info = informations_.at(i);
        const auto index_begin = info.except_indices.cbegin();
        const auto index_end   = info.except_indices.cend();
        const auto chain_begin = info.except_chains.cbegin();
        const auto chain_end   = info.except_chains.cend();

        tmp.clear();
        for(std::size_t j=i+1; j<sys.size(); ++j)
        {
            if(std::find(index_begin, index_end, j) != index_end){continue;}

            const std::size_t j_chain = informations_.at(j).chain_idx;
            if(std::find(chain_begin, chain_end, j_chain) != chain_end)
            {continue;}

            if(length_sq(sys.adjust_direction(sys[j].position - ri)) < rc2)
            {
                tmp.push_back(j);
            }
        }
        this->neighbors_.add_list_for(i, tmp);
    }

    this->current_mergin_ = cutoff_ * mergin_;
    return ;
}

template<typename traitsT>
void VerletList<traitsT>::update(const system_type& sys)
{
    this->current_mergin_ -= sys.max_speed() * dt_ * 2.;
    if(this->current_mergin_ < 0.)
    {
        this->make(sys);
    }
    return ;
}


template<typename traitsT>
void VerletList<traitsT>::update(const system_type& sys, const real_type dt)
{
    this->dt_ = dt;
    this->update(sys);
    return ;
}


} // mjolnir
#endif/* MJOLNIR_CORE_VERLET_LIST */
