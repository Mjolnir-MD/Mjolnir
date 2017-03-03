#ifndef MJOLNIR_UNLIMITED_GRID_CELL_LIST
#define MJOLNIR_UNLIMITED_GRID_CELL_LIST
#include "BoundaryCondition.hpp"
#include "SpatialPartition.hpp"
#include <mjolnir/util/logger.hpp>
#include <functional>
#include <unordered_map>
#include <tuple>
#include <cmath>

namespace std
{

template<>
struct hash<std::tuple<int, int, int>>
{
    typedef std::tuple<int, int, int> argument_type;
    std::size_t operator()(argument_type const& val) const
    {
        return hash<int>()(std::get<0>(val)) ^
               hash<int>()(std::get<1>(val)) ^
               hash<int>()(std::get<2>(val));
    }
};

} // std

namespace mjolnir
{

template<typename traitsT, typename boundaryT = UnlimitedBoundary<traitsT>>
class UnlimitedGridCellList : public SpatialPartition<traitsT>
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
    typedef std::tuple<int, int, int> cell_index_type;
    typedef std::vector<index_list> verlet_list_type;
    typedef std::vector<index_list> except_list_type;

  private:

    constexpr static real_type mesh_epsilon = 1e-6;
    typedef std::unordered_map<cell_index_type, index_list> cell_list_type;

  public:

    UnlimitedGridCellList() = default;
    ~UnlimitedGridCellList() override = default;

    UnlimitedGridCellList(const real_type cutoff, const real_type mergin)
        : dt_(0.), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.),
          inv_cell_size_(1. / (cutoff + mergin + mesh_epsilon))
    {}
    UnlimitedGridCellList(const real_type cutoff, const real_type mergin,
                          const time_type dt)
        : dt_(dt), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.),
          inv_cell_size_(1. / (cutoff + mergin + mesh_epsilon))
    {}

    bool valid() const noexcept override
    {
        return current_mergin_ >= 0. || dt_ == 0.;
    }

    void make  (const particle_container_type& pcon) override;
    void update(const particle_container_type& pcon) override;
    void update(const particle_container_type& pcon,
                const time_type dt) override;

    real_type cutoff() const {return this->cutoff_;}
    real_type mergin() const {return this->mergin_;}

    void set_cutoff(const real_type c);
    void set_mergin(const real_type m);

  private:

    cell_index_type index(const position_type& pos) const;
    cell_index_type
    add(const int x, const int y, const int z, const cell_index_type&) const;

  private:

    real_type dt_;
    real_type cutoff_;
    real_type mergin_;
    real_type current_mergin_;
    real_type inv_cell_size_;

    cell_list_type   cell_list_;

    static Logger& logger_;
};

template<typename traitsT, typename boundaryT>
Logger& UnlimitedGridCellList<traitsT, boundaryT>::logger_ =
        LoggerManager<char>::get_logger("UnlimitedGridCellList");

template<typename traitsT, typename boundaryT>
inline void
UnlimitedGridCellList<traitsT, boundaryT>::set_cutoff(const real_type c)
{
    this->cutoff_ = c;
    this->inv_cell_size_ = 1. / (cutoff_ + mergin_ + mesh_epsilon);
    return;
}

template<typename traitsT, typename boundaryT>
inline void
UnlimitedGridCellList<traitsT, boundaryT>::set_mergin(const real_type m)
{
    this->mergin_ = m;
    this->inv_cell_size_ = 1. / (cutoff_ + mergin_ + mesh_epsilon);
    return;
}

template<typename traitsT, typename boundaryT>
void UnlimitedGridCellList<traitsT, boundaryT>::make(
        const particle_container_type& pcon)
{
    MJOLNIR_LOG_DEBUG("UnlimitedGridCellList<traitsT>::make CALLED");

    this->list_.clear();
    this->list_.resize(pcon.size());

    cell_list_.clear();
    MJOLNIR_LOG_DEBUG("cell_list and verlet_list are cleared");

    const real_type r = cutoff_ + mergin_ + mesh_epsilon;
    MJOLNIR_LOG_DEBUG("unit cell size", r);
    MJOLNIR_LOG_DEBUG("inverse of it", inv_cell_size_);

    std::size_t idx = 0;
    for(auto iter = pcon.cbegin(); iter != pcon.cend(); ++iter)
    {
        const std::tuple<int, int, int> cell = index(iter->position);

        MJOLNIR_LOG_DEBUG(
                idx, "-th particle position", iter->position, "cell index",
                std::get<0>(cell), std::get<1>(cell), std::get<2>(cell));

        cell_list_[cell].push_back(idx);

        // adding its index to all the neighboring cells
        // to make the following code simple
        // (but of cource this consumes memory a lot.
        //  its a memory and performance matter)
        cell_list_[add( 1,  0,  0, cell)].push_back(idx);
        cell_list_[add( 0,  1,  0, cell)].push_back(idx);
        cell_list_[add( 0,  0,  1, cell)].push_back(idx);
        cell_list_[add(-1,  0,  0, cell)].push_back(idx);
        cell_list_[add( 0, -1,  0, cell)].push_back(idx);
        cell_list_[add( 0,  0, -1, cell)].push_back(idx);

        cell_list_[add( 1,  1,  0, cell)].push_back(idx);
        cell_list_[add( 0,  1,  1, cell)].push_back(idx);
        cell_list_[add( 1,  0,  1, cell)].push_back(idx);
        cell_list_[add(-1, -1,  0, cell)].push_back(idx);
        cell_list_[add( 0, -1, -1, cell)].push_back(idx);
        cell_list_[add(-1,  0, -1, cell)].push_back(idx);

        cell_list_[add( 1, -1,  0, cell)].push_back(idx);
        cell_list_[add( 0,  1, -1, cell)].push_back(idx);
        cell_list_[add(-1,  0,  1, cell)].push_back(idx);
        cell_list_[add(-1,  1,  0, cell)].push_back(idx);
        cell_list_[add( 0, -1,  1, cell)].push_back(idx);
        cell_list_[add( 1,  0, -1, cell)].push_back(idx);

        cell_list_[add(-1,  1,  1, cell)].push_back(idx);
        cell_list_[add( 1, -1,  1, cell)].push_back(idx);
        cell_list_[add( 1,  1, -1, cell)].push_back(idx);
        cell_list_[add(-1, -1,  1, cell)].push_back(idx);
        cell_list_[add( 1, -1, -1, cell)].push_back(idx);
        cell_list_[add(-1,  1, -1, cell)].push_back(idx);

        cell_list_[add( 1,  1,  1, cell)].push_back(idx);
        cell_list_[add(-1, -1, -1, cell)].push_back(idx);

        ++idx;
    }
    MJOLNIR_LOG_DEBUG("cell list is updated");

    const real_type r_c = cutoff_ + mergin_;
    const real_type r_c2 = r_c * r_c;

    MJOLNIR_LOG_DEBUG("except list size", this->except_.size());

    MJOLNIR_LOG_DEBUG("exception list is not empty.");
    MJOLNIR_LOG_DEBUG("lookup particles and also except list.");
    for(std::size_t i=0; i<pcon.size(); ++i)
    {
        const coordinate_type ri = pcon[i].position;
        const index_list& clist = cell_list_.at(index(ri));

        MJOLNIR_LOG_DEBUG("particle position", pcon[i].position);
        MJOLNIR_LOG_DEBUG("making verlet list for index", i);
        MJOLNIR_LOG_DEBUG("except list for ", i, "-th value has",
                          this->except_.at(i).size(), "particles");

        const auto cbeg = this->except_.at(i).cbegin();
        const auto cend = this->except_.at(i).cend();
        for(std::size_t j=0; j != clist.size(); ++j)
        {
            MJOLNIR_LOG_DEBUG("looking", j, "-th particle in the cell.",
                              "its index is", clist.at(j));

            std::size_t k = clist.at(j);
            if(k <= i || std::find(cbeg, cend, clist.at(j)) != cend)
                continue;

            if(length_sq(boundary_type::adjust_direction(pcon.at(k).position - ri))
                    < r_c2)
            {
                MJOLNIR_LOG_DEBUG("add index", k, "to verlet list");
                this->list_.at(i).push_back(k);
            }
        }
    }
    this->current_mergin_ = mergin_;

    MJOLNIR_LOG_DEBUG("UnlimitedGridCellList::make() RETURNED");
    return ;
}

template<typename traitsT, typename boundaryT>
void UnlimitedGridCellList<traitsT, boundaryT>::update(
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
inline void UnlimitedGridCellList<traitsT, boundaryT>::update(
        const particle_container_type& pcon, const time_type dt)
{
    this->dt_ = dt;
    this->update(pcon);
    return ;
}

template<typename traitsT, typename boundaryT>
inline typename UnlimitedGridCellList<traitsT, boundaryT>::cell_index_type
UnlimitedGridCellList<traitsT, boundaryT>::index(const position_type& pos) const
{
    return std::make_tuple(static_cast<int>(std::floor(pos[0]*inv_cell_size_)),
                           static_cast<int>(std::floor(pos[1]*inv_cell_size_)),
                           static_cast<int>(std::floor(pos[2]*inv_cell_size_)));
}

template<typename traitsT, typename boundaryT>
inline typename UnlimitedGridCellList<traitsT, boundaryT>::cell_index_type
UnlimitedGridCellList<traitsT, boundaryT>::add(
        const int x, const int y, const int z, const cell_index_type& idx) const
{
    return std::make_tuple(std::get<0>(idx) + x, std::get<1>(idx) + y,
                           std::get<2>(idx) + z);
}

} // mjolnir
#endif/* MJOLNIR_UNLIMITED_GRID_CELL_LIST */
