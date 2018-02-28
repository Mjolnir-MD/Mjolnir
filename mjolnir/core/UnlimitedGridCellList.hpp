#ifndef MJOLNIR_UNLIMITED_GRID_CELL_LIST
#define MJOLNIR_UNLIMITED_GRID_CELL_LIST
#include "BoundaryCondition.hpp"
#include <mjolnir/util/logger.hpp>
#include <functional>
#include <algorithm>
#include <limits>
#include <array>
#include <cmath>

namespace mjolnir
{

template<typename traitsT, std::size_t dimI = 16>
class UnlimitedGridCellList
{
  public:

    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef std::vector<std::size_t> index_array;
    typedef std::vector<index_array> partners_type;

    constexpr static std::size_t dim_size   = dimI;
    constexpr static int         dim        = dim_size;
    constexpr static std::size_t total_size = dim_size * dim_size * dim_size;
    constexpr static real_type mesh_epsilon = 1e-6;
    typedef std::array<int, 3>                        cell_index_type;
    typedef std::array<std::size_t, 26>               neighbor_cell_idx;
    typedef std::pair<index_array, neighbor_cell_idx> unit_cell_type;
    typedef std::array<unit_cell_type, total_size>    cell_list_type;

    struct information
    {
        information() : chain_idx(std::numeric_limits<std::size_t>::max()){}
        std::size_t chain_idx;
        index_array except_chains;
        index_array except_indices;
    };
    typedef std::vector<information> particle_info_type;

  public:

    UnlimitedGridCellList() = default;
    ~UnlimitedGridCellList() = default;
    UnlimitedGridCellList(UnlimitedGridCellList const&) = default;
    UnlimitedGridCellList(UnlimitedGridCellList &&)     = default;
    UnlimitedGridCellList& operator=(UnlimitedGridCellList const&) = default;
    UnlimitedGridCellList& operator=(UnlimitedGridCellList &&)     = default;

    UnlimitedGridCellList(const real_type cutoff, const real_type mergin)
        : dt_(0.), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.),
          inv_cell_size_(1. / (cutoff * (1. + mergin) + mesh_epsilon))
    {}
    UnlimitedGridCellList(const real_type cutoff, const real_type mergin,
                          const real_type dt)
        : dt_(dt), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.),
          inv_cell_size_(1. / (cutoff * (1. + mergin) + mesh_epsilon))
    {}

    bool valid() const noexcept
    {
        return current_mergin_ >= 0. || dt_ == 0.;
    }

    void initialize(const system_type& sys);
    void make  (const system_type& sys);
    void update(const system_type& sys);
    void update(const system_type& sys, const real_type dt);

    real_type cutoff() const {return this->cutoff_;}
    real_type mergin() const {return this->mergin_;}

    void set_cutoff(const real_type c);
    void set_mergin(const real_type m);

    std::size_t& chain_index   (std::size_t i);
    index_array& except_indices(std::size_t i);
    index_array& except_chains (std::size_t i);

    index_array const& partners(std::size_t i) const noexcept {return partners_[i];}

  private:


    std::size_t index(const coordinate_type& pos) const;
    std::size_t index(const cell_index_type& pos) const;
    cell_index_type add(const int x, const int y, const int z,
                        const cell_index_type&) const;

  private:

    real_type dt_;
    real_type cutoff_;
    real_type mergin_;
    real_type current_mergin_;
    real_type inv_cell_size_;
    static Logger& logger_;

    partners_type      partners_;
    particle_info_type informations_;
    cell_list_type     cell_list_;
};

template<typename traitsT, std::size_t N>
Logger& UnlimitedGridCellList<traitsT, N>::logger_ =
        LoggerManager<char>::get_logger("UnlimitedGridCellList");

template<typename traitsT, std::size_t N>
std::size_t& UnlimitedGridCellList<traitsT, N>::chain_index(std::size_t i)
{
    if(this->informations_.size() <= i)
        this->informations_.resize(i+1);
    return this->informations_.at(i).chain_idx;
}

template<typename traitsT, std::size_t N>
typename UnlimitedGridCellList<traitsT, N>::index_array&
UnlimitedGridCellList<traitsT, N>::except_indices(std::size_t i)
{
    if(this->informations_.size() <= i)
        this->informations_.resize(i+1);
    return this->informations_.at(i).except_indices;
}

template<typename traitsT, std::size_t N>
typename UnlimitedGridCellList<traitsT, N>::index_array&
UnlimitedGridCellList<traitsT, N>::except_chains(std::size_t i)
{
    if(this->informations_.size() <= i)
        this->informations_.resize(i+1);
    return this->informations_.at(i).except_chains;
}

template<typename traitsT, std::size_t N>
inline void UnlimitedGridCellList<traitsT, N>::set_cutoff(const real_type c)
{
    this->cutoff_ = c;
    this->inv_cell_size_ = 1. / (cutoff_ * (1. + mergin_) + mesh_epsilon);
    return;
}

template<typename traitsT, std::size_t N>
inline void UnlimitedGridCellList<traitsT, N>::set_mergin(const real_type m)
{
    this->mergin_ = m;
    this->inv_cell_size_ = 1. / (cutoff_ * (1. + mergin_) + mesh_epsilon);
    return;
}

template<typename traitsT, std::size_t N>
void UnlimitedGridCellList<traitsT, N>::make(const system_type& sys)
{
    MJOLNIR_LOG_DEBUG("UnlimitedGridCellList<traitsT>::make CALLED");

    this->partners_.resize(sys.size());
    for(auto& partner : this->partners_) partner.clear();
    for(auto& cell : this->cell_list_) cell.first.clear();
    MJOLNIR_LOG_DEBUG("cell_list and verlet_list are cleared");

    if(informations_.size() < sys.size()) informations_.resize(sys.size());

    std::size_t idx = 0;
    for(const auto& particle : sys)
    {
        MJOLNIR_LOG_DEBUG("set", idx, "-th particle at", index(particle.position));
        cell_list_[index(particle.position)].first.push_back(idx);
        ++idx;
    }
    MJOLNIR_LOG_DEBUG("cell list is updated");

    const real_type r_c = cutoff_ * (1. + mergin_);
    const real_type r_c2 = r_c * r_c;

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const coordinate_type ri   = sys[i].position;
        const unit_cell_type& cell = cell_list_.at(index(ri));

        MJOLNIR_LOG_DEBUG("particle position", sys[i].position);
        MJOLNIR_LOG_DEBUG("cell index", index(ri));
        MJOLNIR_LOG_DEBUG("making verlet list for index", i);

        const auto& info       = informations_[i];
        const auto index_begin = info.except_indices.cbegin();
        const auto index_end   = info.except_indices.cend();
        const auto chain_begin = info.except_chains.cbegin();
        const auto chain_end   = info.except_chains.cend();

        // same cells
        for(std::size_t j : cell.first)
        {
            if(j <= i || std::find(index_begin, index_end, j) != index_end)
                continue;

            const std::size_t j_chain = informations_.at(j).chain_idx;
            if(std::find(chain_begin, chain_end, j_chain) != chain_end)
            {
                continue;
            }

            if(length_sq(sys.adjust_direction(sys.at(j).position - ri)) < r_c2)
            {
                MJOLNIR_LOG_DEBUG("add index", j, "to verlet list", i);
                this->partners_[i].push_back(j);
            }
        }

        // neighbor cells
        for(std::size_t cidx : cell.second)
        {
            MJOLNIR_LOG_DEBUG("neighboring cell index", cidx);

            for(std::size_t j : cell_list_[cidx].first)
            {
                MJOLNIR_LOG_DEBUG("looking particle", j);
                if(j <= i || std::find(index_begin, index_end, j) != index_end)
                    continue;

                const std::size_t j_chain = informations_.at(j).chain_idx;
                if(std::find(chain_begin, chain_end, j_chain) != chain_end)
                    continue;

                if(length_sq(sys.adjust_direction(sys.at(j).position - ri)) < r_c2)
                {
                    MJOLNIR_LOG_DEBUG("add index", j, "to verlet list", i);
                    this->partners_[i].push_back(j);
                }
            }
        }
    }

    for(auto& partner : this->partners_)
    {
        std::sort(partner.begin(), partner.end());
    }

    this->current_mergin_ = cutoff_ * mergin_;

    MJOLNIR_LOG_DEBUG("UnlimitedGridCellList::make() RETURNED");
    return ;
}

template<typename traitsT, std::size_t N>
void UnlimitedGridCellList<traitsT, N>::update(const system_type& sys)
{
    this->current_mergin_ -= sys.max_speed() * dt_ * 2.;

    if(this->current_mergin_ < 0.)
    {
        this->make(sys);
    }
    return ;
}

template<typename traitsT, std::size_t N>
inline void UnlimitedGridCellList<traitsT, N>::update(
        const system_type& sys, const real_type dt)
{
    this->dt_ = dt;
    this->update(sys);
    return ;
}

template<typename traitsT, std::size_t N>
inline std::size_t
UnlimitedGridCellList<traitsT, N>::index(const coordinate_type& pos) const
{
    const int x = static_cast<int>(std::floor(pos[0]*inv_cell_size_)) % dim;
    const int y = static_cast<int>(std::floor(pos[1]*inv_cell_size_)) % dim;
    const int z = static_cast<int>(std::floor(pos[2]*inv_cell_size_)) % dim;
    return index(cell_index_type{{(x < 0) ? x + dim : x,
                                  (y < 0) ? y + dim : y,
                                  (z < 0) ? z + dim : z}});
}

template<typename traitsT, std::size_t N>
inline std::size_t
UnlimitedGridCellList<traitsT, N>::index(const cell_index_type& idx) const
{
    return idx[0] + dim_size * idx[1] + dim_size * dim_size * idx[2];
}


template<typename traitsT, std::size_t N>
inline typename UnlimitedGridCellList<traitsT, N>::cell_index_type
UnlimitedGridCellList<traitsT, N>::add(const int x, const int y, const int z,
        const cell_index_type& idx) const
{
    int ret_x = (idx[0] + x);
    int ret_y = (idx[1] + y);
    int ret_z = (idx[2] + z);
    if(ret_x < 0) ret_x += dim; else if(ret_x >= dim) ret_x -= dim;
    if(ret_y < 0) ret_y += dim; else if(ret_y >= dim) ret_y -= dim;
    if(ret_z < 0) ret_z += dim; else if(ret_z >= dim) ret_z -= dim;
    return cell_index_type{{ret_x, ret_y, ret_z}};
}

template<typename traitsT, std::size_t N>
void UnlimitedGridCellList<traitsT, N>::initialize(const system_type& sys)
{
    MJOLNIR_LOG_DEBUG("UnlimitedGridCellList<traitsT>::initialize CALLED");

    for(auto& cell : this->cell_list_)
    {
        cell.first.reserve(20);
        cell.second.fill(std::numeric_limits<std::size_t>::max());
    }

    for(int x = 0; x < dim; ++x)
    for(int y = 0; y < dim; ++y)
    for(int z = 0; z < dim; ++z)
    {
        const cell_index_type idx{{x, y, z}};
        auto& cell = this->cell_list_[index(idx)];

        MJOLNIR_LOG_DEBUG("cell", x, y, z, "index", index(idx));
        cell.second[ 0] = index(add( 1,  0,  0, idx));
        cell.second[ 1] = index(add( 0,  1,  0, idx));
        cell.second[ 2] = index(add( 0,  0,  1, idx));
        cell.second[ 3] = index(add(-1,  0,  0, idx));
        cell.second[ 4] = index(add( 0, -1,  0, idx));
        cell.second[ 5] = index(add( 0,  0, -1, idx));
        cell.second[ 6] = index(add( 1,  1,  0, idx));
        cell.second[ 7] = index(add( 0,  1,  1, idx));
        cell.second[ 8] = index(add( 1,  0,  1, idx));
        cell.second[ 9] = index(add(-1, -1,  0, idx));
        cell.second[10] = index(add( 0, -1, -1, idx));
        cell.second[11] = index(add(-1,  0, -1, idx));
        cell.second[12] = index(add( 1, -1,  0, idx));
        cell.second[13] = index(add( 0,  1, -1, idx));
        cell.second[14] = index(add(-1,  0,  1, idx));
        cell.second[15] = index(add(-1,  1,  0, idx));
        cell.second[16] = index(add( 0, -1,  1, idx));
        cell.second[17] = index(add( 1,  0, -1, idx));
        cell.second[18] = index(add(-1,  1,  1, idx));
        cell.second[19] = index(add( 1, -1,  1, idx));
        cell.second[20] = index(add( 1,  1, -1, idx));
        cell.second[21] = index(add(-1, -1,  1, idx));
        cell.second[22] = index(add( 1, -1, -1, idx));
        cell.second[23] = index(add(-1,  1, -1, idx));
        cell.second[24] = index(add( 1,  1,  1, idx));
        cell.second[25] = index(add(-1, -1, -1, idx));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1,  0,  0, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 0,  1,  0, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 0,  0,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1,  0,  0, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 0, -1,  0, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 0,  0, -1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1,  1,  0, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 0,  1,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1,  0,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1, -1,  0, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 0, -1, -1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1,  0, -1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1, -1,  0, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 0,  1, -1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1,  0,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1,  1,  0, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 0, -1,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1,  0, -1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1,  1,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1, -1,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1,  1, -1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1, -1,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1, -1, -1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1,  1, -1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1,  1,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1, -1, -1, idx)));
    }
    return;
}



} // mjolnir
#endif/* MJOLNIR_UNLIMITED_GRID_CELL_LIST */
