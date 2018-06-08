#ifndef MJOLNIR_UNLIMITED_GRID_CELL_LIST
#define MJOLNIR_UNLIMITED_GRID_CELL_LIST
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/NeighborList.hpp>
#include <mjolnir/core/ExclusionList.hpp>
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

    typedef ExclusionList exclusion_list_type;
    typedef NeighborList  neighbor_list_type;
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

  public:

    UnlimitedGridCellList()
        : margin_(0.5), current_margin_(-1.0), r_cell_size_(-1.0)
    {}

    ~UnlimitedGridCellList() = default;
    UnlimitedGridCellList(UnlimitedGridCellList const&) = default;
    UnlimitedGridCellList(UnlimitedGridCellList &&)     = default;
    UnlimitedGridCellList& operator=(UnlimitedGridCellList const&) = default;
    UnlimitedGridCellList& operator=(UnlimitedGridCellList &&)     = default;

    UnlimitedGridCellList(const real_type margin)
        : margin_(margin), current_margin_(-1.0), r_cell_size_(-1.0)
    {}

    bool valid() const noexcept
    {
        return current_margin_ >= 0.;
    }

    template<typename PotentialT>
    void initialize(const system_type& sys, const PotentialT& pot);

    template<typename PotentialT>
    void reconstruct(const system_type& sys, const PotentialT& pot)
    {
        this->initialize(sys, pot); // do the same thing as `initialize`
        return;
    }

    void make  (const system_type& sys);
    void update(const system_type& sys);

    real_type cutoff() const noexcept {return this->cutoff_;}
    real_type margin() const noexcept {return this->margin_;}

    range_type partners(std::size_t i) const noexcept {return neighbors_[i];}

  private:

    // calc cell index of the position
    std::size_t calc_index(const coordinate_type& pos) const noexcept
    {
        const auto x =
            static_cast<std::int64_t>(std::floor(pos[0]*r_cell_size_)) % dim;
        const auto y =
            static_cast<std::int64_t>(std::floor(pos[1]*r_cell_size_)) % dim;
        const auto z =
            static_cast<std::int64_t>(std::floor(pos[2]*r_cell_size_)) % dim;

        return calc_index((x<0) ? x+dim : x,
                          (y<0) ? y+dim : y,
                          (z<0) ? z+dim : z);
    }

    std::size_t calc_index(const std::size_t x, const std::size_t y,
                           const std::size_t z) const noexcept
    {
        return x + dim_size * y + dim_size * dim_size * z;
    }

    void set_cutoff(const real_type c) noexcept
    {
        this->cutoff_ = c;
        this->r_cell_size_ = 1 / (cutoff_ * (1 + margin_) * (1+mesh_epsilon));
    }
    void set_margin(const real_type m) noexcept
    {
        this->margin_ = m;
        this->r_cell_size_ = 1 / (cutoff_ * (1 + margin_) * (1+mesh_epsilon));
    }

  private:

    real_type cutoff_;
    real_type margin_;
    real_type current_margin_;
    real_type r_cell_size_;

    exclusion_list_type       exclusion_;
    neighbor_list_type        neighbors_;
    cell_list_type            cell_list_;
    cell_index_container_type index_by_cell_;
    // index_by_cell_ has {particle idx, cell idx} and sorted by cell idx
    // first term of cell list contains first and last idx of index_by_cell
};

template<typename traitsT, std::size_t N>
void UnlimitedGridCellList<traitsT, N>::make(const system_type& sys)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE_DEBUG(UnlimitedGridCellList<traitsT>::make(), 0);

    neighbors_.clear();
    index_by_cell_.resize(sys.size());

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        index_by_cell_[i] = std::make_pair(i, calc_index(sys[i].position));
    }
    std::sort(this->index_by_cell_.begin(), this->index_by_cell_.end(),
        [](const particle_cell_idx_pair& lhs, const particle_cell_idx_pair& rhs)
            noexcept -> bool {return lhs.second < rhs.second;});

    { // assign first and last iterator for each cells
        auto iter = index_by_cell_.cbegin();
        for(std::size_t i=0; i<cell_list_.size(); ++i)
        {
            if(iter == index_by_cell_.cend() || i != iter->second)
            {
                cell_list_[i].first = make_range(iter, iter);
                continue;
            }
            const auto first = iter;
            while(iter != index_by_cell_.cend() && i == iter->second)
            {
                ++iter;
            }
            cell_list_[i].first = make_range(first, iter);
        }
    }

    MJOLNIR_LOG_DEBUG("cell list is updated");

    const real_type r_c  = cutoff_ * (1 + margin_);
    const real_type r_c2 = r_c * r_c;
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        MJOLNIR_SCOPE_DEBUG(for(std::size_t i=0; i<sys.size(); ++i), 0);
        const auto& ri   = sys[i].position;
        const auto& cell = cell_list_.at(this->calc_index(ri));

        MJOLNIR_LOG_DEBUG("particle position", sys[i].position);
        MJOLNIR_LOG_DEBUG("cell index",        calc_index(ri));
        MJOLNIR_LOG_DEBUG("making verlet list for index", i);

        std::vector<std::size_t> partner;
        for(std::size_t cidx : cell.second) // for all adjacent cells...
        {
            for(auto pici : cell_list_[cidx].first)
            {
                const auto j = pici.first;
                MJOLNIR_LOG_DEBUG("looking particle", j);
                if(j <= i || this->exclusion_.is_excluded(i, j))
                {
                    continue;
                }

                const auto& rj = sys[j].position;
                if(length_sq(sys.adjust_direction(rj - ri)) < r_c2)
                {
                    MJOLNIR_LOG_DEBUG("add index", j, "to verlet list", i);
                    partner.push_back(j);
                }
            }
        }
        // make the result consistent with NaivePairCalculation...
        std::sort(partner.begin(), partner.end());
        this->neighbors_.add_list_for(i, partner);
    }

    this->current_margin_ = cutoff_ * margin_;
    return ;
}

template<typename traitsT, std::size_t N>
void UnlimitedGridCellList<traitsT, N>::update(const system_type& sys)
{
    this->current_margin_ -= sys.largest_displacement() * 2;
    if(this->current_margin_ < 0.)
    {
        this->make(sys);
    }
    return ;
}

template<typename traitsT, std::size_t N>
template<typename PotentialT>
void UnlimitedGridCellList<traitsT, N>::initialize(
        const system_type& sys, const PotentialT& pot)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_SCOPE_DEBUG(UnlimitedGridCellList<traitsT>::initialize(), 0);

    this->set_cutoff(pot.max_cutoff_length());
    this->exclusion_.make(sys, pot);

    // initialize cell list
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

    this->make(sys);
    return;
}

} // mjolnir
#endif/* MJOLNIR_UNLIMITED_GRID_CELL_LIST */
