#ifndef MJOLNIR_UNLIMITED_GRID_CELL_LIST
#define MJOLNIR_UNLIMITED_GRID_CELL_LIST
#include <mjolnir/core/NeighborList.hpp>
#include <mjolnir/util/range.hpp>
#include <mjolnir/util/logger.hpp>
#include <functional>
#include <algorithm>
#include <limits>
#include <array>
#include <cmath>

namespace mjolnir
{

template<typename traitsT, std::size_t dimI = 8> // its good to use 2^N as dimI
class UnlimitedGridCellList
{
  public:

    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

    typedef NeighborList           neighbor_list_type;
    typedef neighbor_list_type::range_type range_type;

    constexpr static std::size_t  dim_size  = dimI;
    constexpr static std::int64_t dim       = static_cast<std::int64_t>(dimI);
    constexpr static std::size_t total_size = dim_size * dim_size * dim_size;
    constexpr static real_type mesh_epsilon = 1e-6;

    typedef std::pair<std::size_t, std::size_t> particle_cell_idx_pair;
    typedef std::vector<particle_cell_idx_pair> cell_index_container_type;

    typedef std::array<std::size_t, 27> adjacent_cell_idx;
    typedef std::pair<range<typename cell_index_container_type::const_iterator>,
                      adjacent_cell_idx> cell_type;
    typedef std::array<cell_type, total_size> cell_list_type;

    struct information
    {
        information() : chain_idx(std::numeric_limits<std::size_t>::max()){}
        std::size_t chain_idx;
        std::vector<std::size_t> except_chains;
        std::vector<std::size_t> except_indices;
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
          r_cell_size_(1. / (cutoff * (1. + mergin) + mesh_epsilon))
    {}
    UnlimitedGridCellList(const real_type cutoff, const real_type mergin,
                          const real_type dt)
        : dt_(dt), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.),
          r_cell_size_(1. / (cutoff * (1. + mergin) + mesh_epsilon))
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

    std::size_t&              chain_index   (std::size_t i);
    std::vector<std::size_t>& except_indices(std::size_t i);
    std::vector<std::size_t>& except_chains (std::size_t i);

    range_type partners(std::size_t i) const noexcept {return neighbors_[i];}

  private:

    std::size_t calc_index(const coordinate_type& pos) const noexcept;
    std::size_t calc_index(const std::size_t i, const std::size_t j,
                           const std::size_t k) const noexcept;

  private:

    real_type dt_;
    real_type cutoff_;
    real_type mergin_;
    real_type current_mergin_;
    real_type r_cell_size_;
    static Logger& logger_;

    particle_info_type informations_;
    neighbor_list_type neighbors_;
    cell_list_type     cell_list_;
    cell_index_container_type index_by_cell_;
    // index_by_cell_ has {particle idx, cell idx} and sorted by cell idx
    // first term of cell list contains first and last idx of index_by_cell
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
std::vector<std::size_t>&
UnlimitedGridCellList<traitsT, N>::except_indices(std::size_t i)
{
    if(this->informations_.size() <= i)
        this->informations_.resize(i+1);
    return this->informations_.at(i).except_indices;
}

template<typename traitsT, std::size_t N>
std::vector<std::size_t>&
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
    this->r_cell_size_ = 1. / (cutoff_ * (1. + mergin_) + mesh_epsilon);
    return;
}

template<typename traitsT, std::size_t N>
inline void UnlimitedGridCellList<traitsT, N>::set_mergin(const real_type m)
{
    this->mergin_ = m;
    this->r_cell_size_ = 1. / (cutoff_ * (1. + mergin_) + mesh_epsilon);
    return;
}

template<typename traitsT, std::size_t N>
void UnlimitedGridCellList<traitsT, N>::make(const system_type& sys)
{
    MJOLNIR_LOG_DEBUG("UnlimitedGridCellList<traitsT>::make CALLED");

    neighbors_.clear();
    index_by_cell_.clear();

    if(informations_.size() < sys.size())
    {
        informations_.resize(sys.size());
    }

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        MJOLNIR_LOG_DEBUG(i, "-th particle in", calc_index(sys[i].position));
        index_by_cell_.push_back(
                std::make_pair(i, calc_index(sys[i].position)));
    }
    std::sort(this->index_by_cell_.begin(), this->index_by_cell_.end(),
        [](const std::pair<std::size_t, std::size_t>& lhs,
           const std::pair<std::size_t, std::size_t>& rhs) noexcept -> bool
        {
            return lhs.second < rhs.second;
        });

    { // assign first and last iterator for each cells
        auto iter = index_by_cell_.cbegin();
        for(std::size_t i=0; i<cell_list_.size(); ++i)
        {
            if(i != iter->second)
            {
                cell_list_[i].first = make_range(iter, iter);
                continue;
            }
            const auto first = iter;
            while(i == iter->second)
            {
                ++iter;
            }
            cell_list_[i].first = make_range(first, iter);
        }
    }

    MJOLNIR_LOG_DEBUG("cell list is updated");

    std::vector<std::size_t> tmp;
    const real_type r_c  = cutoff_ * (1. + mergin_);
    const real_type r_c2 = r_c * r_c;
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto& ri   = sys[i].position;
        const auto& cell = cell_list_.at(this->calc_index(ri));

        MJOLNIR_LOG_DEBUG("particle position", sys[i].position);
        MJOLNIR_LOG_DEBUG("cell index", index(ri));
        MJOLNIR_LOG_DEBUG("making verlet list for index", i);

        const auto& info       = informations_[i];
        const auto index_begin = info.except_indices.cbegin();
        const auto index_end   = info.except_indices.cend();
        const auto chain_begin = info.except_chains.cbegin();
        const auto chain_end   = info.except_chains.cend();

        tmp.clear();
        for(std::size_t cidx : cell.second) // for all adjacent cells...
        {
            for(auto pici : cell_list_[cidx].first)
            {
                const auto j = pici.first;
                MJOLNIR_LOG_DEBUG("looking particle", j);
                if(j <= i || std::find(index_begin, index_end, j) != index_end)
                {
                    continue;
                }

                const std::size_t j_chain = informations_.at(j).chain_idx;
                if(std::find(chain_begin, chain_end, j_chain) != chain_end)
                {
                    continue;
                }

                if(length_sq(sys.adjust_direction(sys.at(j).position - ri)) < r_c2)
                {
                    MJOLNIR_LOG_DEBUG("add index", j, "to verlet list", i);
                    tmp.push_back(j);
                }
            }
        }
        // make the result consistent with NaivePairCalculation...
        std::sort(tmp.begin(), tmp.end());
        this->neighbors_.add_list_for(i, tmp);
    }

    this->current_mergin_ = cutoff_ * mergin_;
    MJOLNIR_LOG_DEBUG("UnlimitedGridCellList::make() RETURNED");
    return ;
}

template<typename traitsT, std::size_t N>
void UnlimitedGridCellList<traitsT, N>::update(const system_type& sys)
{
    if(this->current_mergin_ < 0.)
    {
        this->make(sys);
    }
    this->current_mergin_ -= sys.max_speed() * dt_ * 2.;
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
UnlimitedGridCellList<traitsT, N>::calc_index(
        const coordinate_type& pos) const noexcept
{
    const auto x = static_cast<std::int64_t>(std::floor(pos[0]*r_cell_size_)) % dim;
    const auto y = static_cast<std::int64_t>(std::floor(pos[1]*r_cell_size_)) % dim;
    const auto z = static_cast<std::int64_t>(std::floor(pos[2]*r_cell_size_)) % dim;
    return calc_index((x<0) ? x+dim : x, (y<0) ? y+dim : y, (z<0) ? z+dim : z);
}

template<typename traitsT, std::size_t N>
inline std::size_t
UnlimitedGridCellList<traitsT, N>::calc_index(const std::size_t x,
        const std::size_t y, std::size_t z) const noexcept
{
    return x + dim_size * y + dim_size * dim_size * z;
}

template<typename traitsT, std::size_t N>
void UnlimitedGridCellList<traitsT, N>::initialize(const system_type& sys)
{
    MJOLNIR_LOG_DEBUG("UnlimitedGridCellList<traitsT>::initialize CALLED");

    for(int x = 0; x < dim; ++x)
    for(int y = 0; y < dim; ++y)
    for(int z = 0; z < dim; ++z)
    {
        auto& cell = this->cell_list_[calc_index(x, y, z)];

        const std::size_t x_prev = (x ==     0) ? dim-1 : x-1;
        const std::size_t x_next = (x == dim-1) ?     0 : x+1;
        const std::size_t y_prev = (y ==     0) ? dim-1 : y-1;
        const std::size_t y_next = (y == dim-1) ?     0 : y+1;
        const std::size_t z_prev = (z ==     0) ? dim-1 : z-1;
        const std::size_t z_next = (z == dim-1) ?     0 : z+1;

        cell.second[ 0] = calc_index(x_prev, y_prev, z_prev);
        cell.second[ 1] = calc_index(x,      y_prev, z_prev);
        cell.second[ 2] = calc_index(x_next, y_prev, z_prev);
        cell.second[ 3] = calc_index(x_prev, y,      z_prev);
        cell.second[ 4] = calc_index(x,      y,      z_prev);
        cell.second[ 5] = calc_index(x_next, y,      z_prev);
        cell.second[ 6] = calc_index(x_prev, y_next, z_prev);
        cell.second[ 7] = calc_index(x,      y_next, z_prev);
        cell.second[ 8] = calc_index(x_next, y_next, z_prev);

        cell.second[ 9] = calc_index(x_prev, y_prev, z);
        cell.second[10] = calc_index(x,      y_prev, z);
        cell.second[11] = calc_index(x_next, y_prev, z);
        cell.second[12] = calc_index(x_prev, y,      z);
        cell.second[13] = calc_index(x,      y,      z);
        cell.second[14] = calc_index(x_next, y,      z);
        cell.second[15] = calc_index(x_prev, y_next, z);
        cell.second[16] = calc_index(x,      y_next, z);
        cell.second[17] = calc_index(x_next, y_next, z);

        cell.second[18] = calc_index(x_prev, y_prev, z_next);
        cell.second[19] = calc_index(x,      y_prev, z_next);
        cell.second[20] = calc_index(x_next, y_prev, z_next);
        cell.second[21] = calc_index(x_prev, y,      z_next);
        cell.second[22] = calc_index(x,      y,      z_next);
        cell.second[23] = calc_index(x_next, y,      z_next);
        cell.second[24] = calc_index(x_prev, y_next, z_next);
        cell.second[25] = calc_index(x,      y_next, z_next);
        cell.second[26] = calc_index(x_next, y_next, z_next);
    }
    return;
}

} // mjolnir
#endif/* MJOLNIR_UNLIMITED_GRID_CELL_LIST */
