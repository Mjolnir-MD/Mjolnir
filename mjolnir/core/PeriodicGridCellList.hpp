#ifndef MJOLNIR_PERIODIC_GRID_CELL_LIST
#define MJOLNIR_PERIODIC_GRID_CELL_LIST
#include <mjolnir/util/logger.hpp>
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
    typedef std::vector<std::size_t> index_array;
    typedef std::vector<index_array> partners_type;

    constexpr static real_type mesh_epsilon = 1e-6;
    typedef std::array<int, 3>          cell_index_type;
    typedef std::array<std::size_t, 26> neighbor_cell_idx;
    typedef std::pair<index_array, neighbor_cell_idx> unit_cell_type;
    typedef std::vector<unit_cell_type> cell_list_type;

    struct information
    {
        information() : chain_idx(std::numeric_limits<std::size_t>::max()){}
        std::size_t chain_idx;
        index_array except_chains;
        index_array except_indices;
    };
    typedef std::vector<information> particle_info_type;

  public:

    PeriodicGridCellList() = default;
    ~PeriodicGridCellList() = default;
    PeriodicGridCellList(PeriodicGridCellList const&) = default;
    PeriodicGridCellList(PeriodicGridCellList &&)     = default;
    PeriodicGridCellList& operator=(PeriodicGridCellList const&) = default;
    PeriodicGridCellList& operator=(PeriodicGridCellList &&)     = default;

    PeriodicGridCellList(const real_type cutoff, const real_type mergin)
        : dt_(0.), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.),
          inv_cell_size_(1. / (cutoff * (1+mergin) + mesh_epsilon))
    {}

    PeriodicGridCellList(const real_type cutoff, const real_type mergin,
                         const real_type dt)
        : dt_(dt), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.),
          inv_cell_size_(1. / (cutoff * (1+mergin) + mesh_epsilon))
    {}

    bool valid() const noexcept
    {
        return current_mergin_ >= 0. || dt_ == 0.;
    }

    void initialize(const system_type&);
    void make  (const system_type& sys);
    void update(const system_type& sys);
    void update(const system_type& sys, const real_type dt);

    real_type cutoff() const {return this->cutoff_;}
    real_type mergin() const {return this->mergin_;}

    void set_cutoff(const real_type c);
    void set_mergin(const real_type m);

    std::size_t& chain_index   (std::size_t i);
    index_array& except_indices(std::size_t i);
    index_array& except_chains (std::size_t i);

    index_array const& partners(std::size_t i) const noexcept {return partners_[i];}

  private:

    std::size_t index(cell_index_type) const;
    std::size_t index(const coordinate_type& pos) const;
    cell_index_type add(const int x, const int y, const int z,
                        const cell_index_type&) const;

  private:

    real_type   dt_;
    real_type   cutoff_;
    real_type   mergin_;
    real_type   current_mergin_;
    real_type   inv_cell_size_;
    std::size_t dim_x_;
    std::size_t dim_y_;
    std::size_t dim_z_;
    static Logger& logger_;

    coordinate_type    lower_bound_;
    partners_type      partners_;
    particle_info_type informations_;
    cell_list_type     cell_list_;
};

template<typename traitsT>
Logger& PeriodicGridCellList<traitsT>::logger_ =
        LoggerManager<char>::get_logger("PeriodicGridCellList");

template<typename traitsT>
std::size_t& PeriodicGridCellList<traitsT>::chain_index(std::size_t i)
{
    if(this->informations_.size() <= i)
        this->informations_.resize(i+1);
    return this->informations_.at(i).chain_idx;
}

template<typename traitsT>
typename PeriodicGridCellList<traitsT>::index_array&
PeriodicGridCellList<traitsT>::except_indices(std::size_t i)
{
    if(this->informations_.size() <= i)
        this->informations_.resize(i+1);
    return this->informations_.at(i).except_indices;
}

template<typename traitsT>
typename PeriodicGridCellList<traitsT>::index_array&
PeriodicGridCellList<traitsT>::except_chains(std::size_t i)
{
    if(this->informations_.size() <= i)
        this->informations_.resize(i+1);
    return this->informations_.at(i).except_chains;
}

template<typename traitsT>
void PeriodicGridCellList<traitsT>::make(const system_type& sys)
{
    this->lower_bound_ = sys.boundary().lower_bound();
    MJOLNIR_LOG_DEBUG("PeriodicGridCellList<traitsT>::make CALLED");

    this->partners_.resize(sys.size());
    for(auto& partner : this->partners_) partner.clear();
    for(auto& cell : this->cell_list_) cell.first.clear();

    if(informations_.size() < sys.size()) informations_.resize(sys.size());

    std::size_t idx = 0;
    for(const auto& particle : sys)
    {
        MJOLNIR_LOG_DEBUG("set", idx, "-th particle at", index(particle.position));
        cell_list_.at(index(particle.position)).first.push_back(idx);
        ++idx;
    }
    MJOLNIR_LOG_DEBUG("cell list is updated");

    const real_type r_c  = cutoff_ * (1. + mergin_);
    const real_type r_c2 = r_c * r_c;

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const coordinate_type ri = sys[i].position;
        const auto& cell = cell_list_.at(index(ri));

        MJOLNIR_LOG_DEBUG("particle position", sys[i].position);
        MJOLNIR_LOG_DEBUG("making verlet list for index", i);
        MJOLNIR_LOG_DEBUG("except list for ", i, "-th value");

        const auto& info       = informations_[i];
        const auto index_begin = info.except_indices.cbegin();
        const auto index_end   = info.except_indices.cend();
        const auto chain_begin = info.except_chains.cbegin();
        const auto chain_end   = info.except_chains.cend();

        for(std::size_t j : cell.first)
        {
            if(j <= i || std::find(index_begin, index_end, j) != index_end)
                continue;

            const std::size_t j_chain = informations_.at(j).chain_idx;
            if(std::find(chain_begin, chain_end, j_chain) != chain_end)
                continue;

            if(length_sq(sys.adjust_direction(sys.at(j).position - ri)) < r_c2)
            {
                MJOLNIR_LOG_DEBUG("add index", j, "to verlet list of", i);
                this->partners_[i].push_back(j);
            }
        }

        // neighbor cells
        for(std::size_t cidx : cell.second)
        {
            MJOLNIR_LOG_DEBUG("neighboring cell index", cidx);

            for(std::size_t j : cell_list_[cidx].first)
            {
                if(j <= i || std::find(index_begin, index_end, j) != index_end)
                    continue;

                const std::size_t j_chain = informations_.at(j).chain_idx;
                if(std::find(chain_begin, chain_end, j_chain) != chain_end)
                    continue;

                if(length_sq(sys.adjust_direction(sys.at(j).position - ri)) < r_c2)
                {
                    MJOLNIR_LOG_DEBUG("add index", j, "to verlet list of", i);
                    this->partners_[i].push_back(j);
                }
            }
        }
    }

    for(auto& partner : this->partners_)
    {
        std::sort(partner.begin(), partner.end());
    }

    this->current_mergin_ = cutoff_ * mergin_;

    MJOLNIR_LOG_DEBUG("PeriodicGridCellList::make() RETURNED");
    return ;
}



template<typename traitsT>
inline void PeriodicGridCellList<traitsT>::set_cutoff(const real_type c)
{
    this->cutoff_ = c;
    this->inv_cell_size_ = 1. / (cutoff_ * (1. + mergin_) + mesh_epsilon);
    return;
}

template<typename traitsT>
inline void PeriodicGridCellList<traitsT>::set_mergin(const real_type m)
{
    this->mergin_ = m;
    this->inv_cell_size_ = 1. / (cutoff_ * (1. + mergin_) + mesh_epsilon);
    return;
}

template<typename traitsT>
void PeriodicGridCellList<traitsT>::update(const system_type& sys)
{
    this->current_mergin_ -= sys.max_speed() * dt_ * 2.;
    if(this->current_mergin_ < 0.)
        this->make(sys);

    return ;
}

template<typename traitsT>
inline void
PeriodicGridCellList<traitsT>::update(const system_type& sys, const real_type dt)
{
    this->dt_ = dt;
    this->update(sys);
    return ;
}

template<typename traitsT>
inline std::size_t
PeriodicGridCellList<traitsT>::index(const coordinate_type& pos) const
{
    return index(std::array<int, 3>{{
        static_cast<int>(std::floor((pos[0]-lower_bound_[0])*inv_cell_size_)),
        static_cast<int>(std::floor((pos[1]-lower_bound_[1])*inv_cell_size_)),
        static_cast<int>(std::floor((pos[2]-lower_bound_[2])*inv_cell_size_))}});
}

template<typename traitsT>
inline std::size_t
PeriodicGridCellList<traitsT>::index(cell_index_type idx) const
{
    return idx[0] + this->dim_x_ * idx[1] + this->dim_x_ * this->dim_y_ * idx[2];
}


template<typename traitsT>
inline typename PeriodicGridCellList<traitsT>::cell_index_type
PeriodicGridCellList<traitsT>::add(
        const int x, const int y, const int z, const cell_index_type& idx) const
{
    int ret_x = (idx[0] + x);
    int ret_y = (idx[1] + y);
    int ret_z = (idx[2] + z);
    if(ret_x < 0) ret_x += dim_x_; else if(ret_x >= dim_x_) ret_x -= dim_x_;
    if(ret_y < 0) ret_y += dim_y_; else if(ret_y >= dim_y_) ret_y -= dim_y_;
    if(ret_z < 0) ret_z += dim_z_; else if(ret_z >= dim_z_) ret_z -= dim_z_;
    return cell_index_type{{ret_x, ret_y, ret_z}};
}

template<typename traitsT>
void PeriodicGridCellList<traitsT>::initialize(const system_type& sys)
{
    MJOLNIR_LOG_DEBUG("PeriodicGridCellList<traitsT>::initialize CALLED");
    this->lower_bound_ = sys.boundary().lower_bound();
    const auto system_size = sys.boundary().range();
    this->dim_x_ = std::max<std::size_t>(3, std::floor(system_size[0] * inv_cell_size_));
    this->dim_y_ = std::max<std::size_t>(3, std::floor(system_size[1] * inv_cell_size_));
    this->dim_z_ = std::max<std::size_t>(3, std::floor(system_size[2] * inv_cell_size_));
    this->cell_list_.resize(dim_x_ * dim_y_ * dim_z_);

    for(int x = 0; x < dim_x_; ++x)
    for(int y = 0; y < dim_y_; ++y)
    for(int z = 0; z < dim_z_; ++z)
    {
        const cell_index_type idx{{x, y, z}};
        auto& cell = this->cell_list_[index(idx)];

        cell.second[ 0] = index(add( 1,  0,  0, idx));
        cell.second[ 1] = index(add( 0,  1,  0, idx));
        cell.second[ 2] = index(add( 0,  0,  1, idx));
        cell.second[ 3] = index(add(-1,  0,  0, idx));
        cell.second[ 4] = index(add( 0, -1,  0, idx));
        cell.second[ 5] = index(add( 0,  0, -1, idx));
        cell.second[ 6] = index(add( 1,  1,  0, idx));
        cell.second[ 7] = index(add( 0,  1,  1, idx));
        cell.second[ 8] = index(add( 1,  0,  1, idx));
        cell.second[ 9] = index(add(-1, -1,  0, idx));
        cell.second[10] = index(add( 0, -1, -1, idx));
        cell.second[11] = index(add(-1,  0, -1, idx));
        cell.second[12] = index(add( 1, -1,  0, idx));
        cell.second[13] = index(add( 0,  1, -1, idx));
        cell.second[14] = index(add(-1,  0,  1, idx));
        cell.second[15] = index(add(-1,  1,  0, idx));
        cell.second[16] = index(add( 0, -1,  1, idx));
        cell.second[17] = index(add( 1,  0, -1, idx));
        cell.second[18] = index(add(-1,  1,  1, idx));
        cell.second[19] = index(add( 1, -1,  1, idx));
        cell.second[20] = index(add( 1,  1, -1, idx));
        cell.second[21] = index(add(-1, -1,  1, idx));
        cell.second[22] = index(add( 1, -1, -1, idx));
        cell.second[23] = index(add(-1,  1, -1, idx));
        cell.second[24] = index(add( 1,  1,  1, idx));
        cell.second[25] = index(add(-1, -1, -1, idx));

        auto self = std::find(cell.second.cbegin(), cell.second.cend(), index(idx));
        assert(self == cell.second.end());
        auto uniq = std::unique(cell.second.begin(), cell.second.end());
        assert(uniq == cell.second.end());
        for(auto i : cell.second)
        {
            assert(0 <= i && i <= cell_list_.size());
        }
    }
    MJOLNIR_LOG_DEBUG("PeriodicGridCellList<traitsT>::initialize end");
    return;
}


} // mjolnir
#endif /* MJOLNIR_PERIODIC_GRID_CELL_LIST */
