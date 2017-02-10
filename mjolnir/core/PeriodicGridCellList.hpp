#ifndef MJOLNIR_PERIODIC_GRID_CELL_LIST
#define MJOLNIR_PERIODIC_GRID_CELL_LIST
#include "BoundaryCondition.hpp"
#include "SpatialPartition.hpp"
#include <mjolnir/util/logger.hpp>
#include <vector>
#include <cmath>

namespace mjolnir
{

template<typename traitsT, typename boundaryT = PeriodicBoundaryXYZ<traitsT>>
class PeriodicGridCellList : public SpatialPartition<traitsT>
{
  public:

    typedef boundaryT boundary_type;
    typedef SpatialPartition<traitsT> base_type;
    typedef typename base_type::time_type time_type;
    typedef typename base_type::real_type real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::coordinate_type position_type;
    typedef typename base_type::particle_container_type particle_container_type;
    typedef typename base_type::index_type index_type;
    typedef typename base_type::index_list index_list;
    typedef std::array<int, 3>          cell_index_type;
    typedef std::array<index_list*, 26> neighbor_cell_ptr;
    typedef std::pair<index_list, neighbor_cell_ptr> unit_cell_type;
    typedef std::vector<index_list> verlet_list_type;
    typedef std::vector<index_list> except_list_type;
    typedef std::vector<unit_cell_type> cell_list_type;

  private:
    constexpr static real_type mesh_epsilon = 1e-6;

  public:

    PeriodicGridCellList() = default;
    ~PeriodicGridCellList() override = default;

    PeriodicGridCellList(const real_type cutoff, const real_type mergin)
        : dt_(0.), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.),
          inv_cell_size_(1. / (cutoff + mergin + mesh_epsilon))
    {
        initialize();
    }
    PeriodicGridCellList(const real_type cutoff, const real_type mergin,
                         const time_type dt)
        : dt_(dt), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.),
          inv_cell_size_(1. / (cutoff + mergin + mesh_epsilon))
    {
        initialize();
    }

    bool valid() const noexcept override
    {
        return current_mergin_ >= 0. || dt_ == 0.;
    }

    void initialize();
    void make  (const particle_container_type& pcon) override;
    void update(const particle_container_type& pcon) override;
    void update(const particle_container_type& pcon,
                const time_type dt) override;

    real_type const& cutoff() const {return this->cutoff_;}
    real_type const& mergin() const {return this->mergin_;}

    void set_cutoff(const real_type c);
    void set_mergin(const real_type m);

  private:

    std::size_t index(cell_index_type) const;
    std::size_t index(const position_type& pos) const;
    cell_index_type add(const int x, const int y, const int z,
                        const cell_index_type&) const;

  private:

    real_type dt_;
    real_type cutoff_;
    real_type mergin_;
    real_type current_mergin_;
    real_type inv_cell_size_;
    std::size_t dim_x_;
    std::size_t dim_y_;
    std::size_t dim_z_;

    cell_list_type cell_list_;

    static Logger& logger_;
};

template<typename traitsT, typename boundaryT>
Logger& PeriodicGridCellList<traitsT, boundaryT>::logger_ =
        LoggerManager<char>::get_logger("PeriodicGridCellList");

template<typename traitsT, typename boundaryT>
void PeriodicGridCellList<traitsT, boundaryT>::make(
        const particle_container_type& pcon)
{
    MJOLNIR_LOG_DEBUG("PeriodicGridCellList<traitsT>::make CALLED");

    this->list_.clear();
    this->list_.resize(pcon.size());

    for(auto iter = cell_list_.begin(); iter != cell_list_.end(); ++iter)
        iter->first.clear(); // DON'T clear iter->second

    MJOLNIR_LOG_DEBUG("cell_list and verlet_list are cleared");

    std::size_t idx = 0;
    for(auto iter = pcon.cbegin(); iter != pcon.cend(); ++iter)
    {
        cell_list_[index(iter->position - boundary_type::lower_bound())].first.push_back(idx);
        ++idx;
    }
    MJOLNIR_LOG_DEBUG("cell list is updated");

    const real_type r_c  = cutoff_ + mergin_;
    const real_type r_c2 = r_c * r_c;

    MJOLNIR_LOG_DEBUG("except list size", this->except_.size());

    MJOLNIR_LOG_DEBUG("exception list is not empty.");
    MJOLNIR_LOG_DEBUG("lookup particles and also except list.");
    for(std::size_t i=0; i<pcon.size(); ++i)
    {
        const coordinate_type ri = pcon[i].position;
        const auto& cell = cell_list_.at(index(ri - boundary_type::lower_bound()));

        MJOLNIR_LOG_DEBUG("particle position", pcon[i].position);
        MJOLNIR_LOG_DEBUG("making verlet list for index", i);
        MJOLNIR_LOG_DEBUG("except list for ", i, "-th value has",
                          this->except_.at(i).size(), "particles");

        const auto cbeg = this->except_.at(i).cbegin();
        const auto cend = this->except_.at(i).cend();

        for(std::size_t j=0; j<cell.first.size(); ++j)
        {
            MJOLNIR_LOG_DEBUG("looking", j, "-th particle in the cell.",
                              "its index is", clist.at(j));
            const std::size_t k = cell.first.at(j);
            if(k <= i || std::find(cbeg, cend, k) != cend)
                continue;

            if(length_sq(boundary_type::adjust_direction(
                            pcon.at(k).position - ri)) < r_c2)
            {
                MJOLNIR_LOG_DEBUG("add index", k, "to verlet list");
                this->list_.at(i).push_back(k);
            }
        }

        // see neighbor cells
        for(auto iter = cell.second.cbegin(); iter != cell.second.cend(); ++iter)
        {
            for(std::size_t j=0; j<(*iter)->size(); ++j)
            {
                MJOLNIR_LOG_DEBUG("looking", j, "-th particle in the cell.",
                                  "its index is", clist.at(j));
                const std::size_t k = (*iter)->at(j);
                if(k <= i || std::find(cbeg, cend, k) != cend)
                    continue;

                if(length_sq(boundary_type::adjust_direction(
                                pcon.at(k).position - ri)) < r_c2)
                {
                    MJOLNIR_LOG_DEBUG("add index", k, "to verlet list");
                    this->list_.at(i).push_back(k);
                }
            }
        }
    }
    this->current_mergin_ = mergin_;

    MJOLNIR_LOG_DEBUG("PeriodicGridCellList::make() RETURNED");
    return ;
}



template<typename traitsT, typename boundaryT>
inline void
PeriodicGridCellList<traitsT, boundaryT>::set_cutoff(const real_type c)
{
    this->cutoff_ = c;
    this->inv_cell_size_ = 1. / (cutoff_ + mergin_ + mesh_epsilon);
    return;
}

template<typename traitsT, typename boundaryT>
inline void
PeriodicGridCellList<traitsT, boundaryT>::set_mergin(const real_type m)
{
    this->mergin_ = m;
    this->inv_cell_size_ = 1. / (cutoff_ + mergin_ + mesh_epsilon);
    return;
}

template<typename traitsT, typename boundaryT>
void PeriodicGridCellList<traitsT, boundaryT>::update(
        const particle_container_type& pcon)
{
    real_type max_speed = 0.;
    for(auto iter = pcon.cbegin(); iter != pcon.cend(); ++iter)
        max_speed = std::max(max_speed, length_sq(iter->velocity));

    this->current_mergin_ -= std::sqrt(max_speed) * dt_ * 2.;
    if(this->current_mergin_ < 0.)
        this->make(pcon);

    return ;
}

template<typename traitsT, typename boundaryT>
inline void PeriodicGridCellList<traitsT, boundaryT>::update(
        const particle_container_type& pcon, const time_type dt)
{
    this->dt_ = dt;
    this->update(pcon);
    return ;
}

template<typename traitsT, typename boundaryT>
inline std::size_t
PeriodicGridCellList<traitsT, boundaryT>::index(const position_type& pos) const
{
    return index(std::array<int, 3>{{
        static_cast<int>(std::floor(pos[0]*inv_cell_size_)),
        static_cast<int>(std::floor(pos[1]*inv_cell_size_)),
        static_cast<int>(std::floor(pos[2]*inv_cell_size_))}});
}

template<typename traitsT, typename boundaryT>
inline std::size_t
PeriodicGridCellList<traitsT, boundaryT>::index(cell_index_type idx) const
{
    return idx[0] + dim_x_ * idx[1] + dim_x_ * dim_y_ * idx[2];
}


template<typename traitsT, typename boundaryT>
inline typename PeriodicGridCellList<traitsT, boundaryT>::cell_index_type
PeriodicGridCellList<traitsT, boundaryT>::add(
        const int x, const int y, const int z, const cell_index_type& idx) const
{
    int ret_x = (idx[0] + x);
    int ret_y = (idx[1] + y);
    int ret_z = (idx[2] + z);
    if(ret_x < 0) ret_x += dim_x_; else if(ret_x >= dim_x_) ret_x -= dim_x_;
    if(ret_y < 0) ret_y += dim_y_; else if(ret_y >= dim_y_) ret_y -= dim_y_;
    if(ret_z < 0) ret_z += dim_z_; else if(ret_z >= dim_z_) ret_z -= dim_z_;
    return std::array<int, 3>{{ret_x, ret_y, ret_z}};
}

template<typename traitsT, typename boundaryT>
void PeriodicGridCellList<traitsT, boundaryT>::initialize()
{
    this->dim_x_ = boundary_type::system_size()[0] * inv_cell_size_ + 1;
    this->dim_y_ = boundary_type::system_size()[1] * inv_cell_size_ + 1;
    this->dim_z_ = boundary_type::system_size()[2] * inv_cell_size_ + 1;
    this->cell_list_.resize(dim_x_ * dim_y_ * dim_z_);

    for(int x = 0; x < dim_x_; ++x)
    for(int y = 0; y < dim_y_; ++y)
    for(int z = 0; z < dim_z_; ++z)
    {
        const cell_index_type idx{{x, y, z}};
        auto& cell = this->cell_list_[index(idx)];
        cell.second[ 0] = &(cell_list_.at(index(add( 1,  0,  0, idx))).first);
        cell.second[ 1] = &(cell_list_.at(index(add( 0,  1,  0, idx))).first);
        cell.second[ 2] = &(cell_list_.at(index(add( 0,  0,  1, idx))).first);
        cell.second[ 3] = &(cell_list_.at(index(add(-1,  0,  0, idx))).first);
        cell.second[ 4] = &(cell_list_.at(index(add( 0, -1,  0, idx))).first);
        cell.second[ 5] = &(cell_list_.at(index(add( 0,  0, -1, idx))).first);
        cell.second[ 6] = &(cell_list_.at(index(add( 1,  1,  0, idx))).first);
        cell.second[ 7] = &(cell_list_.at(index(add( 0,  1,  1, idx))).first);
        cell.second[ 8] = &(cell_list_.at(index(add( 1,  0,  1, idx))).first);
        cell.second[ 9] = &(cell_list_.at(index(add(-1, -1,  0, idx))).first);
        cell.second[10] = &(cell_list_.at(index(add( 0, -1, -1, idx))).first);
        cell.second[11] = &(cell_list_.at(index(add(-1,  0, -1, idx))).first);
        cell.second[12] = &(cell_list_.at(index(add( 1, -1,  0, idx))).first);
        cell.second[13] = &(cell_list_.at(index(add( 0,  1, -1, idx))).first);
        cell.second[14] = &(cell_list_.at(index(add(-1,  0,  1, idx))).first);
        cell.second[15] = &(cell_list_.at(index(add(-1,  1,  0, idx))).first);
        cell.second[16] = &(cell_list_.at(index(add( 0, -1,  1, idx))).first);
        cell.second[17] = &(cell_list_.at(index(add( 1,  0, -1, idx))).first);
        cell.second[18] = &(cell_list_.at(index(add(-1,  1,  1, idx))).first);
        cell.second[19] = &(cell_list_.at(index(add( 1, -1,  1, idx))).first);
        cell.second[20] = &(cell_list_.at(index(add( 1,  1, -1, idx))).first);
        cell.second[21] = &(cell_list_.at(index(add(-1, -1,  1, idx))).first);
        cell.second[22] = &(cell_list_.at(index(add( 1, -1, -1, idx))).first);
        cell.second[23] = &(cell_list_.at(index(add(-1,  1, -1, idx))).first);
        cell.second[24] = &(cell_list_.at(index(add( 1,  1,  1, idx))).first);
        cell.second[25] = &(cell_list_.at(index(add(-1, -1, -1, idx))).first);
    }
    return;
}


} // mjolnir
#endif /* MJOLNIR_PERIODIC_GRID_CELL_LIST */
