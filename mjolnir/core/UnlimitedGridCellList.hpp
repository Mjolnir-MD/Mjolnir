#ifndef MJOLNIR_CORE_UNLIMITED_GRID_CELL_LIST_HPP
#define MJOLNIR_CORE_UNLIMITED_GRID_CELL_LIST_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/SpatialPartitionBase.hpp>
#include <mjolnir/util/logger.hpp>
#include <functional>
#include <algorithm>
#include <limits>
#include <array>
#include <cmath>

namespace mjolnir
{

template<typename traitsT, typename PotentialT>
class UnlimitedGridCellList final : public SpatialPartitionBase<traitsT, PotentialT>
{
  public:

    using traits_type        = traitsT;
    using potential_type     = PotentialT;
    using base_type          = SpatialPartitionBase<traits_type, potential_type>;

    using system_type        = typename base_type::system_type;
    using boundary_type      = typename base_type::boundary_type;
    using real_type          = typename base_type::real_type;
    using coordinate_type    = typename base_type::coordinate_type;
    using neighbor_list_type = typename base_type::neighbor_list_type;
    using neighbor_type      = typename base_type::neighbor_type;
    using range_type         = typename base_type::range_type;

    static constexpr std::size_t  dim_size()     {return 8u;}
    static constexpr std::int64_t dim()          {return 8; }
    static constexpr std::size_t  total_size()   {return 512u;}
    static constexpr real_type    mesh_epsilon() {return 1e-6;}

    using particle_cell_idx_pair    = std::pair<std::size_t, std::size_t>;
    using cell_index_container_type = std::vector<particle_cell_idx_pair>;
    using cell_index_const_iterator = typename cell_index_container_type::const_iterator;
    using adjacent_cell_idx         = std::array<std::size_t, 27>;
    using cell_type                 = std::pair<range<cell_index_const_iterator>, adjacent_cell_idx>;
    using cell_list_type            = std::array<cell_type, dim() * dim() * dim()>;

  public:

    UnlimitedGridCellList()
        : margin_(0.5), current_margin_(-1.0), r_cell_size_(-1.0)
    {}

    ~UnlimitedGridCellList() override {}
    UnlimitedGridCellList(UnlimitedGridCellList const&) = default;
    UnlimitedGridCellList(UnlimitedGridCellList &&)     = default;
    UnlimitedGridCellList& operator=(UnlimitedGridCellList const&) = default;
    UnlimitedGridCellList& operator=(UnlimitedGridCellList &&)     = default;

    explicit UnlimitedGridCellList(const real_type margin)
        : margin_(margin), current_margin_(-1.0), r_cell_size_(-1.0)
    {}

    bool valid() const noexcept override
    {
        return current_margin_ >= 0.;
    }

    void initialize(neighbor_list_type& neighbors,
                    const system_type& sys, const potential_type& pot) override;
    void make(neighbor_list_type& neighbors,
              const system_type& sys, const potential_type& pot) override;

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
        if(this->current_margin_ < 0)
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
        return new UnlimitedGridCellList(this->margin_);
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
        constexpr std::size_t ds = dim_size();
        return x + ds * y + ds * ds * z;
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

    cell_list_type            cell_list_;
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
void UnlimitedGridCellList<traitsT, potentialT>::make(neighbor_list_type& neighbors,
        const system_type& sys, const potential_type& pot)
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();

    // `participants` is a list that contains indices of particles that are
    // related to the potential.
    const auto& participants = pot.participants();

    neighbors.clear();
    index_by_cell_.resize(participants.size());

    for(std::size_t i=0; i<participants.size(); ++i)
    {
        const auto idx = participants[i];
        index_by_cell_[i] =
            std::make_pair(idx, this->calc_index(sys.position(idx)));
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
                MJOLNIR_LOG_DEBUG("cell ", i, " does not have a particle inside");
                cell_list_[i].first = make_range(iter, iter);
                continue;
            }
            const auto first = iter;
            while(iter != index_by_cell_.cend() && i == iter->second)
            {
                ++iter;
            }
            MJOLNIR_LOG_DEBUG("cell ", i, " has ", std::distance(first, iter),
                              " particle inside");
            cell_list_[i].first = make_range(first, iter);
        }
    }

    MJOLNIR_LOG_DEBUG("cell list is updated");

    const real_type r_c  = cutoff_ * (1 + margin_);
    const real_type r_c2 = r_c * r_c;

    const auto leading_participants = pot.leading_participants();

    std::vector<neighbor_type> partner;
    for(std::size_t idx=0; idx<leading_participants.size(); ++idx)
    {
        partner.clear();
        const auto   i = leading_participants[idx];
        const auto& ri = sys.position(i);

        const auto& cell = cell_list_[this->calc_index(ri)];

        MJOLNIR_LOG_DEBUG("particle position ", sys.position(i));
        MJOLNIR_LOG_DEBUG("cell index ",        calc_index(ri));
        MJOLNIR_LOG_DEBUG("making verlet list for index ", i);

        for(std::size_t cidx : cell.second) // for all adjacent cells...
        {
            MJOLNIR_LOG_DEBUG("neighbor cell index ", cidx);
            for(auto pici : cell_list_[cidx].first)
            {
                const auto j = pici.first;
                MJOLNIR_LOG_DEBUG("looking particle ", j);
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
                    MJOLNIR_LOG_DEBUG("add index ", j, " to list ", i);
                    partner.emplace_back(j, pot.prepare_params(i, j));
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

template<typename traitsT, typename potentialT>
void UnlimitedGridCellList<traitsT, potentialT>::initialize(
        neighbor_list_type& neighbors,
        const system_type& sys, const potential_type& pot)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    constexpr std::int64_t d = dim();

    MJOLNIR_LOG_INFO(pot.name(), " cutoff = ", pot.max_cutoff_length());
    MJOLNIR_LOG_INFO("dimension(independent from system size) = ", d, 'x', d, 'x', d);
    this->set_cutoff(pot.max_cutoff_length());

    // initialize cell list
    for(int x = 0; x < d; ++x)
    for(int y = 0; y < d; ++y)
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
    }

    this->make(neighbors, sys, pot);
    return;
}
} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/forcefield/global/DebyeHuckelPotential.hpp>
#include <mjolnir/forcefield/global/ExcludedVolumePotential.hpp>
#include <mjolnir/forcefield/global/LennardJonesPotential.hpp>
#include <mjolnir/forcefield/global/UniformLennardJonesPotential.hpp>

namespace mjolnir
{
extern template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, DebyeHuckelPotential<SimulatorTraits<double, UnlimitedBoundary>>>;
extern template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, DebyeHuckelPotential<SimulatorTraits<float,  UnlimitedBoundary>>>;

extern template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, ExcludedVolumePotential<SimulatorTraits<double, UnlimitedBoundary>>>;
extern template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, ExcludedVolumePotential<SimulatorTraits<float,  UnlimitedBoundary>>>;

extern template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, LennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>>>;
extern template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, LennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>>>;

extern template class UnlimitedGridCellList<SimulatorTraits<double, UnlimitedBoundary>, UniformLennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>>>;
extern template class UnlimitedGridCellList<SimulatorTraits<float,  UnlimitedBoundary>, UniformLennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>>>;
}
#endif // SEPARATE_BUILD

#endif/* MJOLNIR_UNLIMITED_GRID_CELL_LIST */
