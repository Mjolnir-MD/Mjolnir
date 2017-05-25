#ifndef MJOLNIR_CORE_VERLET_LIST
#define MJOLNIR_CORE_VERLET_LIST
#include <vector>
#include <unordered_map>
#include "SpatialPartition.hpp"
#include "BoundaryCondition.hpp"

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

    VerletList() = default;
    ~VerletList() = default;

    VerletList(const real_type cutoff, const real_type mergin)
        : dt_(0.), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.)
    {}
    VerletList(const real_type cutoff, const real_type mergin,
               const time_type dt)
        : dt_(dt), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.)
    {}

    bool valid() const noexcept
    {
        return current_mergin_ >= 0. || dt_ == 0.;
    }

    void make  (const system_type& sys);
    void update(const system_type& sys);
    void update(const system_type& sys, const real_type dt);

    std::size_t& chain_index   (std::size_t i);
    index_array& except_indices(std::size_t i);
    index_array& except_chains (std::size_t i);

    real_type cutoff() const noexcept {return this->cutoff_;}
    real_type mergin() const noexcept {return this->mergin_;}

    void set_cutoff(const real_type c) noexcept {return this->cutoff_ = c;}
    void set_mergin(const real_type m) noexcept {return this->mergin_ = m;}

    index_array const& partners(std::size_t i) const noexcept {return list_[i];}

  private:

    real_type      dt_;
    real_type      cutoff_;
    real_type      mergin_;
    real_type      current_mergin_;
    static Logger& logger_;

    partners_type      partners_;
    particle_info_type informations_;
};

template<typename traitsT, typename boundaryT>
Logger& VerletList<traitsT, boundaryT>::logger_ =
        LoggerManager<char>::get_logger("VerletList");

template<typename traitsT>
std::size_t& VerletList<traitsT>::chain_index(std::size_t i)
{
    if(this->informations_.size() <= i)
        this->informations_.resize(i+1);
    return this->informations_.at(i).chain_idx;
}

template<typename traitsT>
typename VerletList<traitsT>::index_array&
VerletList<traitsT>::except_indices(std::size_t i)
{
    if(this->informations_.size() <= i)
        this->informations_.resize(i+1);
    return this->informations_.at(i).except_indices;
}

template<typename traitsT>
typename VerletList<traitsT>::index_array&
VerletList<traitsT>::except_chains(std::size_t i)
{
    if(this->informations_.size() <= i)
        this->informations_.resize(i+1);
    return this->informations_.at(i).except_chains;
}

template<typename traitsT>
void VerletList<traitsT>::make(const system_type& sys)
{
    this->partners_.resize(sys.size());
    for(auto& partner : this->partners) partner.clear();

    if(informations_.size() < sys.size()) informations_.resize(sys.size());

    const real_type rc = (cutoff_ + mergin_);
    const real_type rc2 = rc * rc;
    for(std::size_t i=0; i<sys.size()-1; ++i)
    {
        const coordinate_type& ri = sys[i].position;
        const auto index_begin = info.except_indices.cbegin();
        const auto index_end   = info.except_indices.cend();
        const auto chain_begin = info.except_chains.cbegin();
        const auto chain_end   = info.except_chains.cend();

        for(std::size_t j=i+1; j<sys.size(); ++j)
        {
            const std::size_t j_chain = informations_.at(j).chain_idx;
            if(std::find(chain_begin, chain_end, j_chain) != chain_end) continue;
            if(std::find(index_begin, index_end, j)       != index_end) continue;

            if(length_sq(sys.adjust_direction(sys[j].position - ri)) < rc2)
                this->partners_[i].push_back(j);
        }
    }

    this->current_mergin_ = mergin_;
    return ;
}

template<typename traitsT>
void VerletList<traitsT>::update(const system_type& sys)
{
    real_type max_speed = 0.;
    for(const auto& particle : sys)
        max_speed = std::max(max_speed, length_sq(particle.velocity));

    this->current_mergin_ -= std::sqrt(max_speed) * dt_ * 2.;
    if(this->current_mergin_ < 0.)
        this->make(sys);
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
