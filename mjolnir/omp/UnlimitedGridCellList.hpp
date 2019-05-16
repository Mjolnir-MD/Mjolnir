#ifndef MJOLNIR_OMP_UNLIMITED_GRID_CELL_LIST_HPP
#define MJOLNIR_OMP_UNLIMITED_GRID_CELL_LIST_HPP
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/omp/sort.hpp>
#include <mjolnir/core/UnlimitedGridCellList.hpp>

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT,
         typename parameterT, std::size_t dimI>
class UnlimitedGridCellList<
    OpenMPSimulatorTraits<realT, boundaryT>, parameterT, dimI>
{
  public:
    using traits_type         = OpenMPSimulatorTraits<realT, boundaryT>;
    using system_type         = System<traits_type>;
    using real_type           = typename traits_type::real_type;
    using coordinate_type     = typename traits_type::coordinate_type;
    using exclusion_list_type = ExclusionList;
    using parameter_type      = parameterT;
    using neighbor_list_type  = NeighborList<parameter_type>;
    using neighbor_type       = typename neighbor_list_type::neighbor_type;
    using range_type          = typename neighbor_list_type::range_type;

    constexpr static std::size_t  dim_size  = dimI;
    constexpr static std::int64_t dim       = static_cast<std::int64_t>(dimI);
    constexpr static std::size_t total_size = dim_size * dim_size * dim_size;
    constexpr static real_type mesh_epsilon = 1e-6;

    using particle_cell_idx_pair    = std::pair<std::size_t, std::size_t>;
    using cell_index_container_type = std::vector<particle_cell_idx_pair>;
    using cell_index_const_iterator = typename cell_index_container_type::const_iterator;
    using adjacent_cell_idx         = std::array<std::size_t, 27>;
    using cell_type                 = std::pair<range<cell_index_const_iterator>, adjacent_cell_idx>;
    using cell_list_type            = std::array<cell_type, total_size>;

  public:

    UnlimitedGridCellList()
        : margin_(0.5), current_margin_(-1.0), r_cell_size_(-1.0)
    {}

    ~UnlimitedGridCellList() = default;
    UnlimitedGridCellList(UnlimitedGridCellList const&) = default;
    UnlimitedGridCellList(UnlimitedGridCellList &&)     = default;
    UnlimitedGridCellList& operator=(UnlimitedGridCellList const&) = default;
    UnlimitedGridCellList& operator=(UnlimitedGridCellList &&)     = default;

    explicit UnlimitedGridCellList(const real_type margin)
        : margin_(margin), current_margin_(-1.0), r_cell_size_(-1.0)
    {}

    bool valid() const noexcept
    {
        return current_margin_ >= 0.;
    }

    //XXX do NOT call this from parallel region.
    template<typename PotentialT>
    void initialize(const system_type& sys, const PotentialT& pot)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        MJOLNIR_LOG_INFO(pot.name(), " cutoff = ", pot.max_cutoff_length());
        MJOLNIR_LOG_INFO("dimension = ", dim, 'x', dim, 'x', dim);
        this->set_cutoff(pot.max_cutoff_length());
        this->exclusion_.make(sys, pot);

        // initialize cell list
#pragma omp parallel for
        for(int x = 0; x < dim; ++x)
        {
            for(int y = 0; y < dim; ++y)
            {
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
            } // for z
            } // for y
        } // for x (in parallel)

        this->make(sys, pot);
        return;
    }

    //XXX do NOT call this from parallel region
    template<typename PotentialT>
    void make  (const system_type& sys, const PotentialT& pot)
    {
        neighbors_.clear();
        if(index_by_cell_    .size() != sys.size() ||
           index_by_cell_buf_.size() != sys.size())
        {
            index_by_cell_    .resize(sys.size());
            index_by_cell_buf_.resize(sys.size());
        }

#pragma omp parallel for
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            // pairof {particle-idx, cell-idx}
            index_by_cell_[i] = std::make_pair(i, calc_index(sys.position(i)));
        }

        omp::sort(this->index_by_cell_, this->index_by_cell_buf_,
                  [](const particle_cell_idx_pair& lhs,
                     const particle_cell_idx_pair& rhs) noexcept -> bool {
                      return lhs.second < rhs.second;
                  });

        // assign first and last iterator for each cells
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

        const real_type r_c  = cutoff_ * (1 + margin_);
        const real_type r_c2 = r_c * r_c;

//XXX ParallelNeighborList consumes quite a lot of memory resources and makes
//XXX both construction and access slower (especially when system has a small
//XXX number of particles). Because of this, after some benchmarking, I found
//XXX that normal NeighborList works good for most of the cases. Some part of
//XXX neighbor-list construction cannot be parallelized, but it becomes still
//XXX faster.

        std::vector<neighbor_type> partner;
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto& ri   = sys.position(i);
            const auto& cell = cell_list_[this->calc_index(ri)];

            partner.clear();
            for(std::size_t cidx : cell.second) // for all adjacent cells...
            {
                for(auto pici : cell_list_[cidx].first)
                {
                    const auto j = pici.first;
                    if(j <= i || this->exclusion_.is_excluded(i, j))
                    {
                        continue;
                    }

                    const auto& rj = sys.position(j);
                    if(math::length_sq(sys.adjust_direction(rj - ri)) < r_c2)
                    {
                        partner.emplace_back(j, pot.prepare_params(i, j));
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

    //XXX do NOT call this from `parallel` region.
    template<typename PotentialT>
    void update(const real_type   dmargin, const system_type& sys,
                const PotentialT& pot)
    {
        this->current_margin_ -= dmargin;

        if(this->current_margin_ < 0.)
        {
            this->make(sys, pot);
        }
        return ;
    }

    real_type cutoff() const noexcept {return this->cutoff_;}
    real_type margin() const noexcept {return this->margin_;}

    range_type partners(std::size_t i) const noexcept {return neighbors_[i];}

  private:

    // calc cell index of the position
    std::size_t calc_index(const coordinate_type& pos) const noexcept
    {
        const auto x = static_cast<std::int64_t>(std::floor(math::X(pos) * r_cell_size_)) % dim;
        const auto y = static_cast<std::int64_t>(std::floor(math::Y(pos) * r_cell_size_)) % dim;
        const auto z = static_cast<std::int64_t>(std::floor(math::Z(pos) * r_cell_size_)) % dim;

        return this->calc_index(
                (x<0) ? x+dim : x, (y<0) ? y+dim : y, (z<0) ? z+dim : z);
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
    cell_index_container_type index_by_cell_buf_; // buffer for sort
    // index_by_cell_ has {particle idx, cell idx} and sorted by cell idx
    // first term of cell list contains first and last idx of index_by_cell
};

template<typename realT, template<typename, typename> class boundaryT,
         typename parameterT, std::size_t dimI>
constexpr std::size_t UnlimitedGridCellList<
    OpenMPSimulatorTraits<realT, boundaryT>, parameterT, dimI>::dim_size;

template<typename realT, template<typename, typename> class boundaryT,
         typename parameterT, std::size_t dimI>
constexpr std::int64_t UnlimitedGridCellList<
    OpenMPSimulatorTraits<realT, boundaryT>, parameterT, dimI>::dim;

template<typename realT, template<typename, typename> class boundaryT,
         typename parameterT, std::size_t dimI>
constexpr std::size_t UnlimitedGridCellList<
    OpenMPSimulatorTraits<realT, boundaryT>, parameterT, dimI>::total_size;

template<typename realT, template<typename, typename> class boundaryT,
         typename parameterT, std::size_t dimI>
constexpr typename UnlimitedGridCellList<
    OpenMPSimulatorTraits<realT, boundaryT>, parameterT, dimI>::real_type
    UnlimitedGridCellList<OpenMPSimulatorTraits<realT, boundaryT>,
                          parameterT, dimI>::mesh_epsilon;

} // mjolnir
#endif/* MJOLNIR_UNLIMITED_GRID_CELL_LIST */
