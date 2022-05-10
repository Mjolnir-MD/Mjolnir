#ifndef MJOLNIR_OMP_UNLIMITED_GRID_CELL_LIST_HPP
#define MJOLNIR_OMP_UNLIMITED_GRID_CELL_LIST_HPP
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/sort.hpp>
#include <mjolnir/core/UnlimitedGridCellList.hpp>

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT, typename potentialT>
class UnlimitedGridCellList<OpenMPSimulatorTraits<realT, boundaryT>, potentialT>
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
    using parameter_list_type = typename base_type::parameter_list_type;

    static constexpr std::size_t  dim_size () {return 8;}
    static constexpr std::int64_t dim      () {return 8;}
    static constexpr real_type mesh_epsilon() {return 1e-6;}

    using particle_cell_idx_pair    = std::pair<std::size_t, std::size_t>;
    using cell_index_container_type = std::vector<particle_cell_idx_pair>;
    using cell_index_const_iterator = typename cell_index_container_type::const_iterator;
    using adjacent_cell_idx         = std::array<std::size_t, 27>;
    using cell_type                 = std::pair<range<cell_index_const_iterator>, adjacent_cell_idx>;
    using cell_list_type            = std::array<cell_type, dim_size() * dim_size() * dim_size()>;

  public:

    UnlimitedGridCellList()
        : margin_(0.5), current_margin_(-1.0), r_cell_size_(-1.0),
          offsets_threads_(omp_get_max_threads()),
          partners_threads_(omp_get_max_threads()),
          neighbors_threads_(omp_get_max_threads()),
          nranges_threads_(omp_get_max_threads())
    {}

    ~UnlimitedGridCellList() override {};
    UnlimitedGridCellList(UnlimitedGridCellList const&) = default;
    UnlimitedGridCellList(UnlimitedGridCellList &&)     = default;
    UnlimitedGridCellList& operator=(UnlimitedGridCellList const&) = default;
    UnlimitedGridCellList& operator=(UnlimitedGridCellList &&)     = default;

    explicit UnlimitedGridCellList(const real_type margin)
        : margin_(margin), current_margin_(-1.0), r_cell_size_(-1.0),
          offsets_threads_(omp_get_max_threads()),
          partners_threads_(omp_get_max_threads()),
          neighbors_threads_(omp_get_max_threads()),
          nranges_threads_(omp_get_max_threads())
    {}

    bool valid() const noexcept override
    {
        return current_margin_ >= 0.;
    }

    //XXX do NOT call this from parallel region.
    void initialize(neighbor_list_type& neighbors,
                    const system_type& sys, const parameter_list_type& params) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        constexpr std::int64_t d = dim();

        MJOLNIR_LOG_INFO(potential_type::name(), " cutoff = ", params.max_cutoff_length());
        MJOLNIR_LOG_INFO("dimension = ", d, 'x', d, 'x', d);

        this->set_cutoff(params.max_cutoff_length());

        // initialize cell list
#pragma omp parallel for
        for(int x = 0; x < d; ++x)
        {
            for(int y = 0; y < d; ++y)
            {
            for(int z = 0; z < d; ++z)
            {
                auto& cell = this->cell_list_[calc_index(x, y, z)];

                const std::size_t x_prev = (x ==   0) ? d-1 : x-1;
                const std::size_t x_next = (x == d-1) ?   0 : x+1;
                const std::size_t y_prev = (y ==   0) ? d-1 : y-1;
                const std::size_t y_next = (y == d-1) ?   0 : y+1;
                const std::size_t z_prev = (z ==   0) ? d-1 : z-1;
                const std::size_t z_next = (z == d-1) ?   0 : z+1;

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

        this->make(neighbors, sys, params);
        return;
    }

    //XXX do NOT call this from parallel region
    void make(neighbor_list_type& neighbor_list,
              const system_type& sys, const parameter_list_type& params) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
        MJOLNIR_LOG_FUNCTION_DEBUG();

        // `participants` is a list that contains indices of particles that are
        // related to the potential.
        const auto& participants = params.participants();

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

        const auto leading_participants = params.leading_participants();

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
                            partners.emplace_back(j, params.prepare_params(i, j));
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

        principal_neighbors.resize(total_neighbors, neighbor_type{
                std::numeric_limits<std::size_t>::max(),
                typename potential_type::parameter_type{}
            });
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

    //XXX do NOT call this from `parallel` region.
    bool reduce_margin(neighbor_list_type& neighbors, const real_type dmargin,
        const system_type& sys, const parameter_list_type& params) override
    {
        this->current_margin_ -= dmargin;

        if(this->current_margin_ < 0.)
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
        if(this->current_margin_ < 0.)
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
        return new UnlimitedGridCellList(margin_);
    }

  private:

    // calc cell index of the position
    std::size_t calc_index(const coordinate_type& pos) const noexcept
    {
        constexpr std::int64_t d = dim();
        const auto x = std::int64_t(std::floor(math::X(pos) * r_cell_size_)) % d;
        const auto y = std::int64_t(std::floor(math::Y(pos) * r_cell_size_)) % d;
        const auto z = std::int64_t(std::floor(math::Z(pos) * r_cell_size_)) % d;

        return this->calc_index((x<0) ? x+d : x, (y<0) ? y+d : y, (z<0) ? z+d : z);
    }

    std::size_t calc_index(const std::size_t x, const std::size_t y,
                           const std::size_t z) const noexcept
    {
        return x + dim_size() * y + dim_size() * dim_size() * z;
    }

    void set_cutoff(const real_type c) noexcept
    {
        this->cutoff_ = c;
        this->r_cell_size_ = 1 / (cutoff_ * (1 + margin_) * (1+mesh_epsilon()));
    }
    void set_margin(const real_type m) noexcept
    {
        this->margin_ = m;
        this->r_cell_size_ = 1 / (cutoff_ * (1 + margin_) * (1+mesh_epsilon()));
    }

  private:

    real_type cutoff_;
    real_type margin_;
    real_type current_margin_;
    real_type r_cell_size_;

    neighbor_list_type        neighbors_;
    cell_list_type            cell_list_;
    cell_index_container_type index_by_cell_;
    cell_index_container_type index_by_cell_buf_; // buffer for sort
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
#include <mjolnir/forcefield/stoichiometric/GlobalStoichiometricInteractionPotential.hpp>

namespace mjolnir
{
extern template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, DebyeHuckelPotential<double>>;
extern template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, DebyeHuckelPotential<float >>;

extern template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, ExcludedVolumePotential<double>>;
extern template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, ExcludedVolumePotential<float >>;

extern template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, LennardJonesPotential<double>>;
extern template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, LennardJonesPotential<float >>;

extern template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, UniformLennardJonesPotential<double>>;
extern template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, UniformLennardJonesPotential<float >>;

extern template class UnlimitedGridCellList<OpenMPSimulatorTraits<double, UnlimitedBoundary>, GlobalStoichiometricInteractionPotential<double>>;
extern template class UnlimitedGridCellList<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, GlobalStoichiometricInteractionPotential<float >>;
}
#endif // MJOLNIR_SEPARATE_BUILD
#endif/* MJOLNIR_UNLIMITED_GRID_CELL_LIST */
