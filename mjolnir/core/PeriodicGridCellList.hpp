#ifndef MJOLNIR_CORE_PERIODIC_GRID_CELL_LIST_HPP
#define MJOLNIR_CORE_PERIODIC_GRID_CELL_LIST_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/SpatialPartitionBase.hpp>
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
template<typename traitsT, typename potentialT>
class PeriodicGridCellList final
    : public SpatialPartitionBase<traitsT, potentialT>
{
  public:
    using traits_type         = traitsT;
    using potential_type      = potentialT;
    using base_type           = SpatialPartitionBase<traits_type, potential_type>;

    using system_type         = typename base_type::system_type;
    using boundary_type       = typename base_type::boundary_type;
    using real_type           = typename base_type::real_type;
    using coordinate_type     = typename base_type::coordinate_type;
    using neighbor_list_type  = typename base_type::neighbor_list_type;
    using neighbor_type       = typename base_type::neighbor_type;
    using range_type          = typename base_type::range_type;
    using parameter_list_type = typename base_type::parameter_list_type;

    static constexpr real_type mesh_epsilon() {return 1e-6;}

    using particle_cell_idx_pair    = std::pair<std::size_t, std::size_t>;
    using cell_index_container_type = std::vector<particle_cell_idx_pair>;
    using cell_index_const_iterator = typename cell_index_container_type::const_iterator;
    using adjacent_cell_idx         = std::array<std::size_t, 27>;
    using cell_type                 = std::pair<range<cell_index_const_iterator>, adjacent_cell_idx>;
    using cell_list_type            = std::vector<cell_type>;

  public:

    PeriodicGridCellList()
        : cutoff_(0), margin_(1), current_margin_(-1),
          r_x_(-1), r_y_(-1), r_z_(-1), dim_x_(0), dim_y_(0), dim_z_(0)
    {}
    ~PeriodicGridCellList() override {}
    PeriodicGridCellList(PeriodicGridCellList const&) = default;
    PeriodicGridCellList(PeriodicGridCellList &&)     = default;
    PeriodicGridCellList& operator=(PeriodicGridCellList const&) = default;
    PeriodicGridCellList& operator=(PeriodicGridCellList &&)     = default;

    explicit PeriodicGridCellList(const real_type margin)
        : cutoff_(0), margin_(margin), current_margin_(-1),
          r_x_(-1), r_y_(-1), r_z_(-1), dim_x_(0), dim_y_(0), dim_z_(0)
    {}

    bool valid() const noexcept override
    {
        return current_margin_ >= 0.0;
    }

    void initialize(neighbor_list_type& neighbors, const system_type& sys,
                    const parameter_list_type& params) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        const real_type max_cutoff = params.max_cutoff_length();
        this->set_cutoff(max_cutoff);

        MJOLNIR_LOG_INFO(potential_type::name(), " cutoff = ", max_cutoff);

        this->lower_bound_ = sys.boundary().lower_bound();
        this->system_size_ = sys.boundary().width();

        const auto dim_x = std::max<std::size_t>(3, std::floor(math::X(system_size_) * r_x_));
        const auto dim_y = std::max<std::size_t>(3, std::floor(math::Y(system_size_) * r_y_));
        const auto dim_z = std::max<std::size_t>(3, std::floor(math::Z(system_size_) * r_z_));

        // allocate list of cells and set connectivity between them
        this->construct_cells(dim_x, dim_y, dim_z);

        // adjust cell width
        this->r_x_ = real_type(1) / (math::X(this->system_size_) / this->dim_x_);
        this->r_y_ = real_type(1) / (math::Y(this->system_size_) / this->dim_y_);
        this->r_z_ = real_type(1) / (math::Z(this->system_size_) / this->dim_z_);

        // construct neighbor list using cells
        this->make(neighbors, sys, params);
        return;
    }
    void make(neighbor_list_type& neighbors, const system_type& sys,
              const parameter_list_type& params) override;

    bool reduce_margin(neighbor_list_type& neighbors, const real_type dmargin,
                       const system_type& sys, const parameter_list_type& params) override
    {
        this->current_margin_ -= dmargin;
        if(this->current_margin_ < 0)
        {
            this->make(neighbors, sys, params);
            return true;
        }
        return false;
    }
    bool scale_margin(neighbor_list_type& neighbors, const real_type scale,
                const system_type& sys, const parameter_list_type& params) override
    {
        this->current_margin_ = (cutoff_ + current_margin_) * scale - cutoff_;
        if(this->current_margin_ < 0)
        {
            this->make(neighbors, sys, params);
            return true;
        }
        return false;
    }

    real_type cutoff() const noexcept override {return this->cutoff_;}
    real_type margin() const noexcept override {return this->margin_;}

    base_type* clone() const override
    {
        return new PeriodicGridCellList(margin_);
    }

  private:

    std::size_t calc_index(const coordinate_type& pos) const noexcept
    {
        const auto ofs = pos - this->lower_bound_;
        return this->calc_index(// to avoid numeric error
            std::min<std::size_t>(std::floor(math::X(ofs) * this->r_x_), dim_x_-1),
            std::min<std::size_t>(std::floor(math::Y(ofs) * this->r_y_), dim_y_-1),
            std::min<std::size_t>(std::floor(math::Z(ofs) * this->r_z_), dim_z_-1));
    }
    std::size_t calc_index(const std::size_t i, const std::size_t j,
                           const std::size_t k) const noexcept
    {
        return i + this->dim_x_ * j + this->dim_x_ * this->dim_y_ * k;
    }

    // allocate list of cells and set connectivity between them.
    // It does not change cell width, but change the connection between cells.
    void construct_cells(
        const std::size_t dim_x, const std::size_t dim_y, const std::size_t dim_z);

    void set_cutoff(const real_type c) noexcept
    {
        constexpr real_type me = mesh_epsilon();
        this->cutoff_ = c;
        this->r_x_ = 1 / (this->cutoff_ * (1 + this->margin_) + me);
        this->r_y_ = 1 / (this->cutoff_ * (1 + this->margin_) + me);
        this->r_z_ = 1 / (this->cutoff_ * (1 + this->margin_) + me);
        return;
    }
    void set_margin(const real_type m) noexcept
    {
        constexpr real_type me = mesh_epsilon();
        this->margin_ = m;
        this->r_x_ = 1.0 / (this->cutoff_ * (1.0 + this->margin_) + me);
        this->r_y_ = 1.0 / (this->cutoff_ * (1.0 + this->margin_) + me);
        this->r_z_ = 1.0 / (this->cutoff_ * (1.0 + this->margin_) + me);
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
    coordinate_type     system_size_; // size of the boundary condition. in NPT, we need to check it
    cell_list_type      cell_list_;
    cell_index_container_type index_by_cell_;
    // index_by_cell_ has {particle idx, cell idx} and sorted by cell idx
    // first term of cell list contains first and last idx of index_by_cell

#ifdef MJOLNIR_WITH_OPENMP
    // OpenMP implementation uses its own specialization to run it in parallel.
    // So this implementation should not be instanciated with the OpenMP traits.
    static_assert(!is_openmp_simulator_traits<traits_type>::value,
                  "this is the default implementation, not for OpenMP");
#endif
};

template<typename traitsT, typename potentialT>
void PeriodicGridCellList<traitsT, potentialT>::make(
        neighbor_list_type& neighbors, const system_type& sys,
        const parameter_list_type& params)
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();

    // first check the system size because the box size might change under NPT.
    // If it does not match exactly, it means we may need to reconstruct the
    // cell list.
    //
    // If the dimension in each direction is the same, then we don't need to
    // re-construct the cell list because the number of cells and their
    // connectivity are kept.
    {
        const auto current_lower_bound = sys.boundary().lower_bound();
        const auto current_system_size = sys.boundary().width();
        if(math::X(lower_bound_) != math::X(current_lower_bound) ||
           math::Y(lower_bound_) != math::Y(current_lower_bound) ||
           math::Z(lower_bound_) != math::Z(current_lower_bound) ||
           math::X(system_size_) != math::X(current_system_size) ||
           math::Y(system_size_) != math::Y(current_system_size) ||
           math::Z(system_size_) != math::Z(current_system_size) )
        {
            this->lower_bound_ = current_lower_bound;
            this->system_size_ = current_system_size;

            // reset `r_x_`s using cutoff length
            this->set_cutoff(params.max_cutoff_length());

            const auto dim_x = std::max<std::size_t>(3, std::floor(math::X(system_size_) * r_x_));
            const auto dim_y = std::max<std::size_t>(3, std::floor(math::Y(system_size_) * r_y_));
            const auto dim_z = std::max<std::size_t>(3, std::floor(math::Z(system_size_) * r_z_));

            MJOLNIR_LOG_DEBUG("cell = ", current_system_size);
            MJOLNIR_LOG_DEBUG("cell = ", dim_x, "x", dim_y, "x", dim_z);

            // if the number of cells changes, we need to reconstruct the
            // cell size.
            if(dim_x != dim_x_ || dim_y != dim_y_ || dim_z != dim_z_)
            {
                // the number of the cells changes and redefine connections.
                // Also, member variable dim_x_, y_, z_ are updated
                this->construct_cells(dim_x, dim_y, dim_z);
            }

            // update this regardless of the number of cells
            this->r_x_ = 1.0 / (math::X(system_size_) / this->dim_x_);
            this->r_y_ = 1.0 / (math::Y(system_size_) / this->dim_y_);
            this->r_z_ = 1.0 / (math::Z(system_size_) / this->dim_z_);

            MJOLNIR_LOG_DEBUG("reciplocal width of cells in x coordinate = ", r_x_);
            MJOLNIR_LOG_DEBUG("reciplocal width of cells in y coordinate = ", r_y_);
            MJOLNIR_LOG_DEBUG("reciplocal width of cells in z coordinate = ", r_z_);
        }
    }

    // `participants` is a list that contains indices of particles that are
    // related to the potential.
    const auto& participants = params.participants();

    neighbors.clear();
    index_by_cell_.resize(participants.size());

    for(std::size_t i=0; i<participants.size(); ++i)
    {
        const auto idx = participants[i];
        index_by_cell_[i] =
            std::make_pair(idx, this->calc_index(sys.position(idx)));
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

    const auto leading_participants = params.leading_participants();

    std::vector<neighbor_type> partner;
    for(std::size_t idx=0; idx<leading_participants.size(); ++idx)
    {
        partner.clear();
        const auto   i = leading_participants[idx];
        const auto& ri = sys.position(i);

        const auto& cell = cell_list_[calc_index(ri)];

        MJOLNIR_LOG_DEBUG("particle position", sys.position(i));
        MJOLNIR_LOG_DEBUG("making verlet list for index", i);
        MJOLNIR_LOG_DEBUG("except list for ", i, "-th value");

        for(std::size_t cidx : cell.second) // for all adjacent cells...
        {
            for(auto pici : cell_list_[cidx].first)
            {
                const auto j = pici.first;
                MJOLNIR_LOG_DEBUG("looking particle", j);
                if(!params.has_interaction(i, j))
                {
                    continue;
                }
                // here we don't need to search `participants` because
                // cell list contains only participants. non-related
                // particles are already filtered.

                const auto& rj = sys.position(j);
                if(math::length_sq(sys.adjust_direction(ri, rj)) < r_c2)
                {
                    MJOLNIR_LOG_DEBUG("add index", j, "to verlet list", i);
                    partner.emplace_back(j, params.prepare_params(i, j));
                }
            }
        }
        // make the result consistent with NaivePairCalculation...
        std::sort(partner.begin(), partner.end());
        neighbors.add_list_for(i, partner.begin(), partner.end());

        // approximate the memory usage and avoid frequent memory allocation.
        // This block is just for efficiency, and does not affect on the result.
        if(idx == 16)
        {
            neighbors.reserve(neighbors.num_neighbors() / 16,
                              leading_participants.size(), sys.size());
        }
    }

    this->current_margin_ = cutoff_ * margin_;
    return ;
}

// allocate list of cells and set connectivity between them
template<typename traitsT, typename potentialT>
void PeriodicGridCellList<traitsT, potentialT>::construct_cells(
    const std::size_t dim_x, const std::size_t dim_y, const std::size_t dim_z)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    this->dim_x_ = dim_x;
    this->dim_y_ = dim_y;
    this->dim_z_ = dim_z;

    MJOLNIR_LOG_INFO("constructing cells: dimension = ",
                     dim_x_, 'x', dim_y_, 'x', dim_z_);

    if(dim_x_ == 3 || dim_y_ == 3 || dim_z_ == 3)
    {
        MJOLNIR_LOG_WARN("system size is too small (",
                dim_x_, 'x', dim_y_, 'x', dim_z_,
                "). This might cause a problem. Please check the cutoff ratio"
                " and the box size!");
    }

    this->cell_list_.resize(dim_x_ * dim_y_ * dim_z_);

    const int dimx = static_cast<int>(dim_x_);
    const int dimy = static_cast<int>(dim_y_);
    const int dimz = static_cast<int>(dim_z_);

    for(int x = 0; x < dimx; ++x)
    {
    for(int y = 0; y < dimy; ++y)
    {
    for(int z = 0; z < dimz; ++z)
    {
        auto& cell = this->cell_list_[calc_index(x, y, z)];

        const std::size_t x_prev = (x ==        0) ? dim_x_ - 1 : x - 1;
        const std::size_t x_next = (x == dimx - 1) ?          0 : x + 1;
        const std::size_t y_prev = (y ==        0) ? dim_y_ - 1 : y - 1;
        const std::size_t y_next = (y == dimy - 1) ?          0 : y + 1;
        const std::size_t z_prev = (z ==        0) ? dim_z_ - 1 : z - 1;
        const std::size_t z_next = (z == dimz - 1) ?          0 : z + 1;

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
            assert(i <= cell_list_.size());
        }
    }
    }
    }
    return;
}

} // mjolnir

// #ifdef MJOLNIR_SEPARATE_BUILD
// #include <mjolnir/forcefield/global/DebyeHuckelPotential.hpp>
// #include <mjolnir/forcefield/global/ExcludedVolumePotential.hpp>
// #include <mjolnir/forcefield/global/LennardJonesPotential.hpp>
//
// namespace mjolnir
// {
// extern template class PeriodicGridCellList<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
// extern template class PeriodicGridCellList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
//
// extern template class PeriodicGridCellList<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
// extern template class PeriodicGridCellList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
//
// extern template class PeriodicGridCellList<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
// extern template class PeriodicGridCellList<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
// }
// #endif // SEPARATE_BUILD

#endif /* MJOLNIR_PERIODIC_GRID_CELL_LIST */
