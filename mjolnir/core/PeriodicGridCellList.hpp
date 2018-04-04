#ifndef MJOLNIR_PERIODIC_GRID_CELL_LIST
#define MJOLNIR_PERIODIC_GRID_CELL_LIST
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/NeighborList.hpp>
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/util/range.hpp>
#include <mjolnir/util/logger.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

namespace mjolnir
{

// XXX: almost same as UnlimitedGridCellList.
// the difference between UnlimitedGridCellList is only the number of cells.
// PeriodicGridCellList can optimize the number of cells using boundary size.
template<typename traitsT>
class PeriodicGridCellList
{
  public:

    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

    typedef ExclusionList exclusion_list_type;
    typedef NeighborList  neighbor_list_type;
    typedef neighbor_list_type::range_type range_type;

    constexpr static real_type mesh_epsilon = 1e-6;

    typedef std::pair<std::size_t, std::size_t> particle_cell_idx_pair;
    typedef std::vector<particle_cell_idx_pair> cell_index_container_type;

    typedef std::array<std::size_t, 27> adjacent_cell_idx;
    typedef std::pair<range<typename cell_index_container_type::const_iterator>,
                      adjacent_cell_idx> cell_type;
    typedef std::vector<cell_type> cell_list_type;

  public:

    PeriodicGridCellList()
        : mergin_(1), current_mergin_(-1), r_x_(-1), r_y_(-1), r_z_(-1)
    {}
    ~PeriodicGridCellList() = default;
    PeriodicGridCellList(PeriodicGridCellList const&) = default;
    PeriodicGridCellList(PeriodicGridCellList &&)     = default;
    PeriodicGridCellList& operator=(PeriodicGridCellList const&) = default;
    PeriodicGridCellList& operator=(PeriodicGridCellList &&)     = default;

    PeriodicGridCellList(const real_type mergin)
        : mergin_(mergin), current_mergin_(-1), r_x_(-1), r_y_(-1), r_z_(-1)
    {}

    bool valid() const noexcept
    {
        return current_mergin_ >= 0.0;
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

    real_type cutoff() const {return this->cutoff_;}
    real_type mergin() const {return this->mergin_;}

    // after calling this, neighbor list should be reconstructed!
    void set_cutoff(const real_type c);
    void set_mergin(const real_type m);

    range_type partners(std::size_t i) const noexcept {return neighbors_[i];}

  private:

    std::size_t calc_index(const coordinate_type& pos) const noexcept
    {
        return calc_index(
            static_cast<std::size_t>(std::floor((pos[0]-lower_bound_[0])*r_x_)),
            static_cast<std::size_t>(std::floor((pos[1]-lower_bound_[1])*r_y_)),
            static_cast<std::size_t>(std::floor((pos[2]-lower_bound_[2])*r_z_)));
    }

    std::size_t calc_index(const std::size_t i, const std::size_t j,
                           const std::size_t k) const noexcept
    {
        return i + this->dim_x_ * j + this->dim_x_ * this->dim_y_ * k;
    }

  private:

    real_type   cutoff_;
    real_type   mergin_;
    real_type   current_mergin_;
    real_type   r_x_;
    real_type   r_y_;
    real_type   r_z_;
    std::size_t dim_x_;
    std::size_t dim_y_;
    std::size_t dim_z_;
    static Logger& logger_;

    coordinate_type     lower_bound_;
    neighbor_list_type  neighbors_;
    exclusion_list_type exclusion_;
    cell_list_type      cell_list_;
    cell_index_container_type index_by_cell_;
    // index_by_cell_ has {particle idx, cell idx} and sorted by cell idx
    // first term of cell list contains first and last idx of index_by_cell
};

template<typename traitsT>
Logger& PeriodicGridCellList<traitsT>::logger_ =
        LoggerManager<char>::get_logger("PeriodicGridCellList");

template<typename traitsT>
void PeriodicGridCellList<traitsT>::make(const system_type& sys)
{
    MJOLNIR_LOG_DEBUG("PeriodicGridCellList<traitsT>::make CALLED");

    neighbors_.clear();
    index_by_cell_.resize(sys.size());

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        index_by_cell_[i] = std::make_pair(i, calc_index(sys[i].position));
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

    const real_type r_c  = cutoff_ * (1. + mergin_);
    const real_type r_c2 = r_c * r_c;
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto& ri = sys[i].position;
        const auto& cell = cell_list_.at(calc_index(ri));

        MJOLNIR_LOG_DEBUG("particle position", sys[i].position);
        MJOLNIR_LOG_DEBUG("making verlet list for index", i);
        MJOLNIR_LOG_DEBUG("except list for ", i, "-th value");

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

                if(length_sq(sys.adjust_direction(sys.at(j).position - ri)) < r_c2)
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

    this->current_mergin_ = cutoff_ * mergin_;
    MJOLNIR_LOG_DEBUG("PeriodicGridCellList::make() RETURNED");
    return ;
}



template<typename traitsT>
inline void PeriodicGridCellList<traitsT>::set_cutoff(const real_type c)
{
    this->cutoff_ = c;
    this->r_x_ = 1.0 / (this->cutoff_ * (1.0 + this->mergin_) + mesh_epsilon);
    this->r_y_ = 1.0 / (this->cutoff_ * (1.0 + this->mergin_) + mesh_epsilon);
    this->r_z_ = 1.0 / (this->cutoff_ * (1.0 + this->mergin_) + mesh_epsilon);
    return;
}

template<typename traitsT>
inline void PeriodicGridCellList<traitsT>::set_mergin(const real_type m)
{
    this->mergin_ = m;
    this->r_x_ = 1.0 / (this->cutoff_ * (1.0 + this->mergin_) + mesh_epsilon);
    this->r_y_ = 1.0 / (this->cutoff_ * (1.0 + this->mergin_) + mesh_epsilon);
    this->r_z_ = 1.0 / (this->cutoff_ * (1.0 + this->mergin_) + mesh_epsilon);
    return;
}

template<typename traitsT>
void PeriodicGridCellList<traitsT>::update(const system_type& sys)
{
    // TODO consider boundary size
    this->current_mergin_ -= sys.largest_displacement() * 2.;
    if(this->current_mergin_ < 0.)
    {
        this->make(sys);
    }
    return ;
}

template<typename traitsT>
template<typename PotentialT>
void PeriodicGridCellList<traitsT>::initialize(
        const system_type& sys, const PotentialT& pot)
{
    MJOLNIR_LOG_DEBUG("PeriodicGridCellList<traitsT>::initialize CALLED");
    this->set_cutoff(pot.max_cutoff_length());
    this->exclusion_.make(sys, pot);

    this->lower_bound_ = sys.boundary().lower_bound();
    const auto system_size = sys.boundary().width();

    this->dim_x_ = std::max<std::size_t>(3, std::floor(system_size[0] * r_x_));
    this->dim_y_ = std::max<std::size_t>(3, std::floor(system_size[1] * r_y_));
    this->dim_z_ = std::max<std::size_t>(3, std::floor(system_size[2] * r_z_));

    if(dim_x_ == 3 || dim_y_ == 3 || dim_z_ == 3)
    {
        std::cerr << "WARNING: cell size might be too small: number of grids =("
                  << dim_x_ << ", " << dim_y_ << ", " << dim_z_ << ")\n";
    }

    // it may expand cell a bit (to fit system range)
    this->r_x_ = 1.0 / (system_size[0] / this->dim_x_);
    this->r_y_ = 1.0 / (system_size[1] / this->dim_y_);
    this->r_z_ = 1.0 / (system_size[2] / this->dim_z_);

    this->cell_list_.resize(dim_x_ * dim_y_ * dim_z_);

    for(int x = 0; x < dim_x_; ++x)
    for(int y = 0; y < dim_y_; ++y)
    for(int z = 0; z < dim_z_; ++z)
    {
        auto& cell = this->cell_list_[calc_index(x, y, z)];

        const std::size_t x_prev = (x ==        0) ? dim_x_-1 : x-1;
        const std::size_t x_next = (x == dim_x_-1) ?        0 : x+1;
        const std::size_t y_prev = (y ==        0) ? dim_y_-1 : y-1;
        const std::size_t y_next = (y == dim_y_-1) ?        0 : y+1;
        const std::size_t z_prev = (z ==        0) ? dim_z_-1 : z-1;
        const std::size_t z_next = (z == dim_z_-1) ?        0 : z+1;

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

        auto uniq = std::unique(cell.second.begin(), cell.second.end());
        assert(uniq == cell.second.end());
        for(auto i : cell.second)
        {
            assert(0 <= i && i <= cell_list_.size());
        }
    }
    this->make(sys);
    MJOLNIR_LOG_DEBUG("PeriodicGridCellList<traitsT>::initialize end");
    return;
}


} // mjolnir
#endif /* MJOLNIR_PERIODIC_GRID_CELL_LIST */
