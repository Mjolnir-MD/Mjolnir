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
template<typename traitsT, typename parameterT>
class PeriodicGridCellList
{
  public:
    using traits_type         = traitsT;
    using system_type         = System<traits_type>;
    using real_type           = typename traits_type::real_type;
    using coordinate_type     = typename traits_type::coordinate_type;
    using exclusion_list_type = ExclusionList;
    using parameter_type      = parameterT;
    using neighbor_list_type  = NeighborList<parameter_type>;
    using neighbor_type       = typename neighbor_list_type::neighbor_type;
    using range_type          = typename neighbor_list_type::range_type;

    constexpr static real_type mesh_epsilon = 1e-6;

    using particle_cell_idx_pair    = std::pair<std::size_t, std::size_t>;
    using cell_index_container_type = std::vector<particle_cell_idx_pair>;
    using cell_index_const_iterator = typename cell_index_container_type::const_iterator;
    using adjacent_cell_idx         = std::array<std::size_t, 27>;
    using cell_type                 = std::pair<range<cell_index_const_iterator>, adjacent_cell_idx>;
    using cell_list_type            = std::vector<cell_type>;

  public:

    PeriodicGridCellList()
        : cutoff_(0), margin_(1), current_margin_(-1),
          r_x_(-1), r_y_(-1), r_z_(-1), dim_x_(0), dim_y_(0), dim_z_(0),
    {}
    ~PeriodicGridCellList() = default;
    PeriodicGridCellList(PeriodicGridCellList const&) = default;
    PeriodicGridCellList(PeriodicGridCellList &&)     = default;
    PeriodicGridCellList& operator=(PeriodicGridCellList const&) = default;
    PeriodicGridCellList& operator=(PeriodicGridCellList &&)     = default;

    explicit PeriodicGridCellList(const real_type margin)
        : cutoff_(0), margin_(margin), current_margin_(-1),
          r_x_(-1), r_y_(-1), r_z_(-1), dim_x_(0), dim_y_(0), dim_z_(0),
    {}

    bool valid() const noexcept
    {
        return current_margin_ >= 0.0;
    }

    template<typename PotentialT>
    void initialize(const system_type& sys, const PotentialT& pot);

    template<typename PotentialT>
    void reconstruct(const system_type& sys, const PotentialT& pot)
    {
        this->initialize(sys, pot); // do the same thing as `initialize`
        return;
    }

    template<typename PotentialT>
    void make  (const system_type& sys, const PotentialT& pot);

    template<typename PotentialT>
    void update(const real_type, const system_type&, const PotentialT&);

    real_type cutoff() const {return this->cutoff_;}
    real_type margin() const {return this->margin_;}

    range_type partners(std::size_t i) const noexcept {return neighbors_[i];}

  private:

    std::size_t calc_index(const coordinate_type& pos) const noexcept
    {
        const auto ofs = pos - this->lower_bound_;
        return this->calc_index(std::floor(math::X(ofs) * this->r_x_),
                                std::floor(math::Y(ofs) * this->r_y_),
                                std::floor(math::Z(ofs) * this->r_z_));
    }

    std::size_t calc_index(const std::size_t i, const std::size_t j,
                           const std::size_t k) const noexcept
    {
        return i + this->dim_x_ * j + this->dim_x_ * this->dim_y_ * k;
    }

    void set_cutoff(const real_type c) noexcept
    {
        this->cutoff_ = c;
        this->r_x_ = 1 / (this->cutoff_ * (1 + this->margin_) + mesh_epsilon);
        this->r_y_ = 1 / (this->cutoff_ * (1 + this->margin_) + mesh_epsilon);
        this->r_z_ = 1 / (this->cutoff_ * (1 + this->margin_) + mesh_epsilon);
        return;
    }
    void set_margin(const real_type m) noexcept
    {
        this->margin_ = m;
        this->r_x_ = 1.0 / (this->cutoff_ * (1.0 + this->margin_) + mesh_epsilon);
        this->r_y_ = 1.0 / (this->cutoff_ * (1.0 + this->margin_) + mesh_epsilon);
        this->r_z_ = 1.0 / (this->cutoff_ * (1.0 + this->margin_) + mesh_epsilon);
        return;
    }

  private:

    real_type   cutoff_;
    real_type   margin_;
    real_type   current_margin_;
    real_type   r_x_;
    real_type   r_y_;
    real_type   r_z_;
    std::size_t dim_x_;
    std::size_t dim_y_;
    std::size_t dim_z_;

    coordinate_type     lower_bound_;
    neighbor_list_type  neighbors_;
    exclusion_list_type exclusion_;
    cell_list_type      cell_list_;
    cell_index_container_type index_by_cell_;
    // index_by_cell_ has {particle idx, cell idx} and sorted by cell idx
    // first term of cell list contains first and last idx of index_by_cell
};

template<typename traitsT, typename parameterT>
template<typename PotentialT>
void PeriodicGridCellList<traitsT, parameterT>::update(
        const real_type dmargin, const system_type& sys, const PotentialT& pot)
{
    this->current_margin_ -= dmargin;
    if(this->current_margin_ < 0.)
    {
        this->make(sys, pot);
    }
    return ;
}

template<typename traitsT, typename parameterT>
template<typename PotentialT>
void PeriodicGridCellList<traitsT, parameterT>::make(
        const system_type& sys, const PotentialT& pot)
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();

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

    const real_type r_c  = cutoff_ * (1. + margin_);
    const real_type r_c2 = r_c * r_c;
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto& ri = sys[i].position;
        const auto& cell = cell_list_[calc_index(ri)];

        MJOLNIR_LOG_DEBUG("particle position", sys[i].position);
        MJOLNIR_LOG_DEBUG("making verlet list for index", i);
        MJOLNIR_LOG_DEBUG("except list for ", i, "-th value");

        std::vector<neighbor_type> partner;
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

                if(math::length_sq(sys.adjust_direction(sys[j].position - ri)) < r_c2)
                {
                    MJOLNIR_LOG_DEBUG("add index", j, "to verlet list", i);
                    partner.emplace_back(j, pot.prepair_params(i, j));
                }
            }
        }
        // make the result consistent with NaivePairCalculation...
        std::sort(partner.begin(), partner.end());
        this->neighbors_.add_list_for(i, partner.begin(), partner.end());
    }

    this->current_margin_ = cutoff_ * margin_;
    return ;
}

template<typename traitsT, typename parameterT>
template<typename PotentialT>
void PeriodicGridCellList<traitsT, parameterT>::initialize(
        const system_type& sys, const PotentialT& pot)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const real_type max_cutoff = pot.max_cutoff_length();
    this->set_cutoff(max_cutoff);
    this->exclusion_.make(sys, pot);

    MJOLNIR_LOG_INFO(pot.name(), " cutoff = ", max_cutoff);

    this->lower_bound_ = sys.boundary().lower_bound();
    const auto system_size = sys.boundary().width();

    this->dim_x_ = std::max<std::size_t>(3, std::floor(math::X(system_size) * r_x_));
    this->dim_y_ = std::max<std::size_t>(3, std::floor(math::Y(system_size) * r_y_));
    this->dim_z_ = std::max<std::size_t>(3, std::floor(math::Z(system_size) * r_z_));

    MJOLNIR_LOG_INFO("dimension = ", dim_x_, 'x', dim_y_, 'x', dim_z_);

    if(dim_x_ == 3 || dim_y_ == 3 || dim_z_ == 3)
    {
        MJOLNIR_LOG_WARN("system size is too small (",
                dim_x_, 'x', dim_y_, 'x', dim_z_,
                "). This might cause a problem. Please check the cutoff ratio"
                " and the box size!");
    }

    // it may expand cell a bit (to fit system range)
    this->r_x_ = 1.0 / (math::X(system_size) / this->dim_x_);
    this->r_y_ = 1.0 / (math::Y(system_size) / this->dim_y_);
    this->r_z_ = 1.0 / (math::Z(system_size) / this->dim_z_);

    MJOLNIR_LOG_DEBUG("reciplocal width of cells in x coordinate = ", r_x_);
    MJOLNIR_LOG_DEBUG("reciplocal width of cells in y coordinate = ", r_y_);
    MJOLNIR_LOG_DEBUG("reciplocal width of cells in z coordinate = ", r_z_);

    this->cell_list_.resize(dim_x_ * dim_y_ * dim_z_);

    for(int x = 0; x < dim_x_; ++x)
    {
    for(int y = 0; y < dim_y_; ++y)
    {
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
    }
    }
    this->make(sys, pot);
    return;
}


} // mjolnir
#endif /* MJOLNIR_PERIODIC_GRID_CELL_LIST */
