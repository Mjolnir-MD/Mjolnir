#ifndef MJOLNIR_AXIS_ALIGNED_PLANE_HPP
#define MJOLNIR_AXIS_ALIGNED_PLANE_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{

/*! @brief axis aligned Plane for ExternalDistanceInteraction.                *
 *  @details represents flat Plane. It provides a method to calculate         *
 *           a distance between particle and the plane, force direction       *
 *           a position of a particle. It also provide a functionality of     *
 *           a neighbor-list. NormalAxis means the index of axis (0 means x,  *
 *           1 means y, 2 means z).                                           */
template<typename traitsT, std::size_t NormalAxis>
class AxisAlignedPlane
{
  public:
    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type        real_type;
    typedef typename traits_type::coordinate_type  coordinate_type;
    typedef typename traits_type::boundary_type    boundary_type;
    constexpr static std::size_t axis = NormalAxis;
    static_assert(axis < 3, "0 <= AxisAlignedPlane::axis < 3");

  public:

    AxisAlignedPlane(const real_type position, const real_type margin = 1)
        : position_(position), margin_(margin), current_margin_(-1)
    {}

    //XXX can be negative! because normal vector can define the direction...
    real_type calc_distance(
            const coordinate_type& pos, const boundary_type& bd) const
    {
        coordinate_type ref(0,0,0);
        ref[axis] = this->position_;

        return bd.adjust_direction(pos - ref)[axis];
    }

    coordinate_type calc_force_direction(
            const coordinate_type& pos, const boundary_type& bd) const
    {
        coordinate_type ref(0,0,0);
        ref[axis] = this->position_;

        const real_type sign = std::copysign(
                real_type(1.0), bd.adjust_direction(pos - ref)[axis]);

        coordinate_type f(0,0,0);
        f[axis] = sign;
        return f;
    }

    template<typename Potential>
    void initialize(const system_type& sys, const Potential& pot)
    {
        // assuming pot is already initialized!
        this->cutoff_         = pot.max_cutoff_length();
        this->participant_    = pot.participants();
        this->current_margin_ = 0.0;
        this->make(sys);
        return;
    }

    //! update exclusion list and cutoff length.
    template<typename Potential>
    void reconstruct(const system_type& sys, const Potential& pot)
    {
        this->initialize(sys, pot);
        return;
    }

    // neighbor-list stuff
    void make  (const system_type& sys);
    void update(const system_type& sys);

    std::vector<std::size_t> const& neighbors() const noexcept
    {return this->neighbors_;}

  private:

    real_type position_; // position in the axis.
    real_type cutoff_, margin_, current_margin_;
    std::vector<std::size_t> neighbors_;   // being inside of cutoff range
    std::vector<std::size_t> participant_; // particle that interacts with
};

template<typename traitsT, std::size_t NormalAxis>
constexpr std::size_t AxisAlignedPlane<traitsT, NormalAxis>::axis;

template<typename traitsT, std::size_t NormalAxis>
void AxisAlignedPlane<traitsT, NormalAxis>::make(const system_type& sys)
{
    this->neighbors_.clear();
    const real_type threshold = this->cutoff_ * (1 + this->margin_);
    const real_type thr2 = threshold * threshold;

    for(std::size_t i : this->participant_)
    {
        const auto d = this->calc_distance(sys[i].position, sys.boundary());
        if(d * d <= thr2)
        {
            this->neighbors_.push_back(i);
        }
    }

    this->current_margin_ = this->cutoff_ * this->margin_;
    return;
}

template<typename traitsT, std::size_t NormalAxis>
void AxisAlignedPlane<traitsT, NormalAxis>::update(const system_type& sys)
{
    this->current_margin_ -= sys.largest_displacement();
    if(this->current_margin_ < 0)
    {
        this->make(sys);
    }
    return;
}

} // mjolnir
#endif // MJOLNIR_AXIS_ALIGNED_PLANE_HPP
