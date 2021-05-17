#ifndef MJOLNIR_OMP_PERIODIC_GRID_CELL_LIST_HPP
#define MJOLNIR_OMP_PERIODIC_GRID_CELL_LIST_HPP
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/sort.hpp>
#include <mjolnir/core/PeriodicGridCellList.hpp>

#include <omp.h>

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT,
         typename potentialT>
class PeriodicGridCellList<OpenMPSimulatorTraits<realT, boundaryT>, potentialT>
    final : public SpatialPartitionBase<OpenMPSimulatorTraits<realT, boundaryT>, potentialT>
{
  public:

    using traits_type        = OpenMPSimulatorTraits<realT, boundaryT>;
    using potential_type     = potentialT;
    using base_type          = SpatialPartitionBase<traits_type, potential_type>;

    using system_type        = typename base_type::system_type;
    using boundary_type      = typename base_type::boundary_type;
    using real_type          = typename base_type::real_type;
    using coordinate_type    = typename base_type::coordinate_type;
    using neighbor_list_type = typename base_type::neighbor_list_type;
    using neighbor_type      = typename base_type::neighbor_type;
    using range_type         = typename base_type::range_type;

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
          r_x_(-1), r_y_(-1), r_z_(-1), dim_x_(0), dim_y_(0), dim_z_(0),
          offsets_threads_(omp_get_max_threads()),
          partners_threads_(omp_get_max_threads()),
          neighbors_threads_(omp_get_max_threads()),
          nranges_threads_(omp_get_max_threads())
    {}
    ~PeriodicGridCellList() override {}
    PeriodicGridCellList(PeriodicGridCellList const&) = default;
    PeriodicGridCellList(PeriodicGridCellList &&)     = default;
    PeriodicGridCellList& operator=(PeriodicGridCellList const&) = default;
    PeriodicGridCellList& operator=(PeriodicGridCellList &&)     = default;

    explicit PeriodicGridCellList(const real_type margin)
        : cutoff_(0), margin_(margin), current_margin_(-1),
          r_x_(-1), r_y_(-1), r_z_(-1), dim_x_(0), dim_y_(0), dim_z_(0),
          offsets_threads_(omp_get_max_threads()),
          partners_threads_(omp_get_max_threads()),
          neighbors_threads_(omp_get_max_threads()),
          nranges_threads_(omp_get_max_threads())
    {}

    bool valid() const noexcept override
    {
        return current_margin_ >= 0.0;
    }

    void initialize(neighbor_list_type& neighbors,
                    const system_type& sys, const potential_type& pot) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        const real_type max_cutoff = pot.max_cutoff_length();
        this->set_cutoff(max_cutoff);

        MJOLNIR_LOG_INFO(pot.name(), " cutoff = ", max_cutoff);

        this->lower_bound_ = sys.boundary().lower_bound();
        this->system_size_ = sys.boundary().width();

        const auto dim_x = std::max<std::size_t>(3, std::floor(math::X(this->system_size_) * r_x_));
        const auto dim_y = std::max<std::size_t>(3, std::floor(math::Y(this->system_size_) * r_y_));
        const auto dim_z = std::max<std::size_t>(3, std::floor(math::Z(this->system_size_) * r_z_));

        // allocate list of cells and set connectivity between them
        this->construct_cells(dim_x, dim_y, dim_z);

        // adjust cell width
        this->r_x_ = real_type(1) / (math::X(this->system_size_) / this->dim_x_);
        this->r_y_ = real_type(1) / (math::Y(this->system_size_) / this->dim_y_);
        this->r_z_ = real_type(1) / (math::Z(this->system_size_) / this->dim_z_);

        // construct neighbor list using cells
        this->make(neighbors, sys, pot);
        return;
    }

    void make(neighbor_list_type& neighbor_list,
              const system_type& sys, const potential_type& pot) override
    {
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
                this->set_cutoff(pot.max_cutoff_length());

                const auto dim_x = std::max<std::size_t>(3, std::floor(math::X(system_size_) * r_x_));
                const auto dim_y = std::max<std::size_t>(3, std::floor(math::Y(system_size_) * r_y_));
                const auto dim_z = std::max<std::size_t>(3, std::floor(math::Z(system_size_) * r_z_));

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
        const auto& participants = pot.participants();

        neighbor_list.clear();
        if(index_by_cell_    .size() != participants.size() ||
           index_by_cell_buf_.size() != participants.size())
        {
            index_by_cell_    .resize(participants.size());
            index_by_cell_buf_.resize(participants.size());
        }

#pragma omp parallel for
        for(std::size_t i=0; i<participants.size(); ++i)
        {
            const auto idx = participants[i];
            index_by_cell_[i] =
                std::make_pair(idx, this->calc_index(sys.position(idx)));
        }

        omp::sort(this->index_by_cell_, this->index_by_cell_buf_,
            [](const particle_cell_idx_pair& lhs,
               const particle_cell_idx_pair& rhs) noexcept -> bool
            {
                return lhs.second < rhs.second;
            });

#pragma omp parallel for
        for(std::size_t cell_idx=0; cell_idx<cell_list_.size(); ++cell_idx)
        {
            auto iter = std::find_if(index_by_cell_.cbegin(), index_by_cell_.cend(),
                [cell_idx](const particle_cell_idx_pair& item) noexcept -> bool {
                    return item.second == cell_idx;
                });
            if(iter == index_by_cell_.cend())
            {
                cell_list_[cell_idx].first = make_range(iter, iter);
                continue;
            }
            const auto first = iter; // the first iter of the range
            while(iter != index_by_cell_.cend() && iter->second == cell_idx)
            {
                ++iter;
            }
            cell_list_[cell_idx].first = make_range(first, iter);
        }

        const real_type r_c  = cutoff_ * (1. + margin_);
        const real_type r_c2 = r_c * r_c;

        const auto leading_participants = pot.leading_participants();

        assert(std::is_sorted(leading_participants.begin(), leading_participants.end()));

        constexpr std::size_t nil = std::numeric_limits<std::size_t>::max();
        std::fill(offsets_threads_.begin(), offsets_threads_.end(), nil);

#pragma omp parallel
        {
            const std::size_t num_threads = omp_get_num_threads();
            const std::size_t thread_id   = omp_get_thread_num();

            const std::size_t dN = leading_participants.size() / num_threads;
            const std::size_t first = dN *  thread_id;
            const std::size_t last  = (thread_id+1 == num_threads) ?
                leading_participants.size() : dN * (thread_id+1);

            const std::size_t first_idx = leading_participants[first];
            const std::size_t  last_idx = (thread_id+1 == num_threads) ?
                leading_participants[leading_participants.size()-1] + 1 : leading_participants[last];

            this->offsets_threads_[thread_id] = first_idx;

            auto& partners  = this->partners_threads_[thread_id];
            auto& neighbors = this->neighbors_threads_[thread_id];
            auto& nranges   = this->nranges_threads_[thread_id];

            neighbors.clear(); // keep capacity
            nranges  .clear();
            nranges.resize(last_idx + 1 - first_idx, 0);

            for(std::size_t idx=first; idx<last; ++idx)
            {
                partners.clear();
                const auto     i = leading_participants[idx];
                const auto    ri = sys.position(i);
                const auto& cell = cell_list_[calc_index(ri)];

                for(std::size_t cidx : cell.second) // for all adjacent cells...
                {
                    for(auto pici : cell_list_[cidx].first)
                    {
                        const auto j = pici.first;
                        if(!pot.has_interaction(i, j))
                        {
                            continue;
                        }

                        // here we don't need to search `participants` because
                        // cell list contains only participants. non-related
                        // particles are already filtered.

                        const auto& rj = sys.position(j);
                        if(math::length_sq(sys.adjust_direction(ri, rj)) < r_c2)
                        {
                            partners.emplace_back(j, pot.prepare_params(i, j));
                        }
                    }
                }
                // make the result consistent with NaivePairCalculation...
                std::sort(partners.begin(), partners.end());

                nranges[i - first_idx    ] = neighbors.size();
                nranges[i - first_idx + 1] = neighbors.size() + partners.size();

                neighbors.reserve(neighbors.size() + partners.size());
                std::copy(partners.begin(), partners.end(),
                          std::back_inserter(neighbors));

                if(idx == first+16)
                {
                    neighbors.reserve(dN * neighbors.size() / 16);
                }
            }
        }

        auto& principal_neighbors = neighbor_list.neighbors();
        auto& principal_ranges    = neighbor_list.ranges();

        std::size_t total_neighbors = 0;
        for(std::size_t th=0; th < offsets_threads_.size(); ++th)
        {
            const auto offset = offsets_threads_[th];
            if(offset == std::numeric_limits<std::size_t>::max())
            {
                break;
            }
            total_neighbors += neighbors_threads_[th].size();
        }

        principal_neighbors.resize(total_neighbors);
        principal_ranges   .resize(sys.size(), 0);

#pragma omp parallel for schedule(static, 1)
        for(std::size_t th=0; th < offsets_threads_.size(); ++th)
        {
            const auto index_offset = offsets_threads_[th];
            if(index_offset != std::numeric_limits<std::size_t>::max())
            {
                std::size_t neighbor_offset = 0;
                for(std::size_t t=0; t<th; ++t)
                {
                    neighbor_offset += neighbors_threads_[t].size();
                }
                auto& neighbors = this->neighbors_threads_[th];
                std::copy(neighbors.begin(), neighbors.end(),
                          principal_neighbors.begin() + neighbor_offset);

                auto& nranges = this->nranges_threads_[th];
                std::transform(nranges.begin() + 1, nranges.end(),
                      principal_ranges.begin() + 1 + index_offset,
                      [=](const std::size_t i) {return i + neighbor_offset;});
            }
        }
        this->current_margin_ = cutoff_ * margin_;
        return ;
    }

    bool reduce_margin(neighbor_list_type& neighbors, const real_type dmargin,
        const system_type& sys, const potential_type& pot) override
    {
        this->current_margin_ -= dmargin;
        if(this->current_margin_ < 0)
        {
            this->make(neighbors, sys, pot);
            return true;
        }
        return false;
    }
    bool scale_margin(neighbor_list_type& neighbors, const real_type scale,
        const system_type& sys, const potential_type& pot) override
    {
        this->current_margin_ = (cutoff_ + current_margin_) * scale - cutoff_;
        if(this->current_margin_ < 0.)
        {
            this->make(neighbors, sys, pot);
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
        return this->calc_index(
            std::min<std::size_t>(std::floor(math::X(ofs) * this->r_x_), dim_x_-1),
            std::min<std::size_t>(std::floor(math::Y(ofs) * this->r_y_), dim_y_-1),
            std::min<std::size_t>(std::floor(math::Z(ofs) * this->r_z_), dim_z_-1));
    }

    std::size_t calc_index(const std::size_t i, const std::size_t j,
                           const std::size_t k) const noexcept
    {
        return i + this->dim_x_ * j + this->dim_x_ * this->dim_y_ * k;
    }

    void set_cutoff(const real_type c) noexcept
    {
        this->cutoff_ = c;
        this->r_x_ = 1 / (this->cutoff_ * (1 + this->margin_) + mesh_epsilon());
        this->r_y_ = 1 / (this->cutoff_ * (1 + this->margin_) + mesh_epsilon());
        this->r_z_ = 1 / (this->cutoff_ * (1 + this->margin_) + mesh_epsilon());
        return;
    }
    void set_margin(const real_type m) noexcept
    {
        this->margin_ = m;
        this->r_x_ = 1.0 / (this->cutoff_ * (1.0 + this->margin_) + mesh_epsilon());
        this->r_y_ = 1.0 / (this->cutoff_ * (1.0 + this->margin_) + mesh_epsilon());
        this->r_z_ = 1.0 / (this->cutoff_ * (1.0 + this->margin_) + mesh_epsilon());
        return;
    }

    void construct_cells(const std::size_t dim_x, const std::size_t dim_y,
                         const std::size_t dim_z)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->dim_x_ = dim_x;
        this->dim_y_ = dim_y;
        this->dim_z_ = dim_z;

        MJOLNIR_LOG_INFO("dimension = ", dim_x_, 'x', dim_y_, 'x', dim_z_);

        if(dim_x_ == 3 || dim_y_ == 3 || dim_z_ == 3)
        {
            MJOLNIR_LOG_WARN("system size is too small (",
                    dim_x_, 'x', dim_y_, 'x', dim_z_,
                    "). This might cause a problem. Please check the cutoff ratio"
                    " and the box size!");
        }

        this->cell_list_.resize(dim_x_ * dim_y_ * dim_z_);

        const int dimx = dim_x_;
        const int dimy = dim_y_;
        const int dimz = dim_z_;

#pragma omp parallel for
        for(int x = 0; x < dimx; ++x)
        {
            for(int y = 0; y < dimy; ++y)
            {
                for(int z = 0; z < dimz; ++z)
                {
                    auto& cell = this->cell_list_[calc_index(x, y, z)];

                    const std::size_t x_prev = (x ==        0) ? dimx - 1 : x - 1;
                    const std::size_t x_next = (x == dimx - 1) ?        0 : x + 1;
                    const std::size_t y_prev = (y ==        0) ? dimy - 1 : y - 1;
                    const std::size_t y_next = (y == dimy - 1) ?        0 : y + 1;
                    const std::size_t z_prev = (z ==        0) ? dimz - 1 : z - 1;
                    const std::size_t z_next = (z == dimz - 1) ?        0 : z + 1;

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

    coordinate_type           lower_bound_;
    coordinate_type           system_size_; // size of the boundary condition. in NPT, we need to check it
    cell_list_type            cell_list_;
    cell_index_container_type index_by_cell_;
    cell_index_container_type index_by_cell_buf_;
    // index_by_cell_ has {particle idx, cell idx} and sorted by cell idx
    // first term of cell list contains first and last idx of index_by_cell

    std::vector<std::size_t> offsets_threads_;
    std::vector<std::vector<neighbor_type>> partners_threads_;
    std::vector<std::vector<neighbor_type>> neighbors_threads_;
    std::vector<std::vector<std::size_t>>   nranges_threads_;
};
} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/forcefield/global/DebyeHuckelPotential.hpp>
#include <mjolnir/forcefield/global/ExcludedVolumePotential.hpp>
#include <mjolnir/forcefield/global/LennardJonesPotential.hpp>
#include <mjolnir/forcefield/global/UniformLennardJonesPotential.hpp>
#include <mjolnir/forcefield/global/HardCoreExcludedVolumePotential.hpp>
#include <mjolnir/forcefield/global/InversePowerPotential.hpp>
#include <mjolnir/forcefield/global/WCAPotential.hpp>

namespace mjolnir
{
extern template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, HardCoreExcludedVolumePotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, HardCoreExcludedVolumePotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, InversePowerPotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, InversePowerPotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class PeriodicGridCellList<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, WCAPotential<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class PeriodicGridCellList<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, WCAPotential<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
}
#endif // SEPARATE_BUILD

#endif /* MJOLNIR_PERIODIC_GRID_CELL_LIST */
