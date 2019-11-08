#ifndef MJOLNIR_OMP_PERIODIC_GRID_CELL_LIST_HPP
#define MJOLNIR_OMP_PERIODIC_GRID_CELL_LIST_HPP
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/sort.hpp>
#include <mjolnir/core/PeriodicGridCellList.hpp>

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

    void initialize(neighbor_list_type& neighbors,
                    const system_type& sys, const potential_type& pot) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        const real_type max_cutoff = pot.max_cutoff_length();
        this->set_cutoff(max_cutoff);

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

        const int dimx = dim_x_;
        const int dimy = dim_y_;
        const int dimz = dim_x_;

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
        this->make(neighbors, sys, pot);
        return;
    }

    void make(neighbor_list_type& neighbors,
              const system_type& sys, const potential_type& pot) override
    {
        // `participants` is a list that contains indices of particles that are
        // related to the potential.
        const auto& participants = pot.participants();

        neighbors.clear();
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

//XXX ParallelNeighborList consumes quite a lot of memory resources and makes
//XXX both construction and access slower (especially when system has a small
//XXX number of particles). Because of this, after some benchmarking, I found
//XXX that normal NeighborList works good for most of the cases. Some part of
//XXX neighbor-list construction cannot be parallelized, but it becomes still
//XXX faster.
        const auto leading_participants = pot.leading_participants();

        std::vector<neighbor_type> partner;
        for(std::size_t idx=0; idx<leading_participants.size(); ++idx)
        {
            partner.clear();
            const auto   i = leading_participants[idx];
            const auto& ri = sys.position(i);

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
                    if(math::length_sq(sys.adjust_direction(rj - ri)) < r_c2)
                    {
                        partner.emplace_back(j, pot.prepare_params(i, j));
                    }
                }
            }
            // make the result consistent with NaivePairCalculation...
            std::sort(partner.begin(), partner.end());
            neighbors.add_list_for(i, partner.begin(), partner.end());
        }

        this->current_margin_ = cutoff_ * margin_;
        return ;
    }


    void update(neighbor_list_type& neighbors, const real_type dmargin,
                const system_type& sys, const potential_type& pot) override
    {
        this->current_margin_ -= dmargin;
        if(this->current_margin_ < 0)
        {
            this->make(neighbors, sys, pot);
        }
        return ;
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
    cell_list_type            cell_list_;
    cell_index_container_type index_by_cell_;
    cell_index_container_type index_by_cell_buf_;
    // index_by_cell_ has {particle idx, cell idx} and sorted by cell idx
    // first term of cell list contains first and last idx of index_by_cell
};
} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/potential/global/DebyeHuckelPotential.hpp>
#include <mjolnir/potential/global/ExcludedVolumePotential.hpp>
#include <mjolnir/potential/global/LennardJonesPotential.hpp>
#include <mjolnir/potential/global/UniformLennardJonesPotential.hpp>

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
}
#endif // SEPARATE_BUILD

#endif /* MJOLNIR_PERIODIC_GRID_CELL_LIST */
