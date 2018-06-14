#ifndef MJOLNIR_AXIS_ALIGNED_PLANE_HPP
#define MJOLNIR_AXIS_ALIGNED_PLANE_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{

/* These classes represents Normal Vectors. (+/- * XYZ) */
template<typename traitsT>
class PositiveXDirection
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;

    static constexpr std::size_t index = 0; // index in coordinate_type
    static constexpr real_type   sign  = 1.0;
    static coordinate_type invoke(const real_type length = 1.0) noexcept
    {return coordinate_type(length, 0, 0);}
};
template<typename T>
constexpr std::size_t PositiveXDirection<T>::index;
template<typename T>
constexpr typename PositiveXDirection<T>::real_type PositiveXDirection<T>::sign;

template<typename traitsT>
class NegativeXDirection
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;

    static constexpr std::size_t index =  0;
    static constexpr real_type   sign  = -1.0;
    static coordinate_type invoke(const real_type length = 1.0) noexcept
    {return coordinate_type(-length, 0, 0);}
};
template<typename T>
constexpr std::size_t NegativeXDirection<T>::index;
template<typename T>
constexpr typename NegativeXDirection<T>::real_type NegativeXDirection<T>::sign;

template<typename traitsT>
class PositiveYDirection
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;

    static constexpr std::size_t index = 1; // index in coordinate_type
    static constexpr real_type   sign  = 1.0;
    static coordinate_type invoke(const real_type length = 1.0) noexcept
    {return coordinate_type(0, length, 0);}
};
template<typename T>
constexpr std::size_t PositiveYDirection<T>::index;
template<typename T>
constexpr typename PositiveYDirection<T>::real_type PositiveYDirection<T>::sign;

template<typename traitsT>
class NegativeYDirection
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;

    static constexpr std::size_t index =  1;
    static constexpr real_type   sign  = -1.0;
    static coordinate_type invoke(const real_type length = 1.0) noexcept
    {return coordinate_type(0, -length, 0);}
};
template<typename T>
constexpr std::size_t NegativeYDirection<T>::index;
template<typename T>
constexpr typename NegativeYDirection<T>::real_type NegativeYDirection<T>::sign;

template<typename traitsT>
class PositiveZDirection
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;

    static constexpr std::size_t index = 2; // index in coordinate_type
    static constexpr real_type   sign  = 1.0;
    static coordinate_type invoke(const real_type length = 1.0) noexcept
    {return coordinate_type(0, 0, length);}
};
template<typename T>
constexpr std::size_t PositiveZDirection<T>::index;
template<typename T>
constexpr typename PositiveZDirection<T>::real_type PositiveZDirection<T>::sign;

template<typename traitsT>
class NegativeZDirection
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;

    static constexpr std::size_t index =  2;
    static constexpr real_type   sign  = -1.0;
    static coordinate_type invoke(const real_type length = 1.0) noexcept
    {return coordinate_type(0, 0, -length);}
};
template<typename T>
constexpr std::size_t NegativeZDirection<T>::index;
template<typename T>
constexpr typename NegativeZDirection<T>::real_type NegativeZDirection<T>::sign;

/*! @brief axis aligned Plane for ExternalDistanceInteraction.                *
 *  @details represents flat Plane. It provides a method to calculate         *
 *           a distance between particle and the plane, force direction       *
 *           a position of a particle. It also provide a functionality of     *
 *           a neighbor-list.                                                 */
template<typename traitsT, template<typename> class normal_axisT>
class AxisAlignedPlane
{
  public:
    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type        real_type;
    typedef typename traits_type::coordinate_type  coordinate_type;
    typedef typename traits_type::boundary_type    boundary_type;
    typedef normal_axisT<traits_type>              normal_axis_type;
    static constexpr std::size_t axis_index = normal_axis_type::index;
    static constexpr real_type   axis_sign  = normal_axis_type::sign;

  public:

    AxisAlignedPlane(const real_type position, const real_type margin = 1)
        : origin_(0, 0, 0), margin_(margin), current_margin_(-1)
    {
        this->origin_[axis_index] = position;
    }

    //XXX can be negative! because normal vector can define the direction...
    inline real_type
    calc_distance(const coordinate_type& pos, const boundary_type& bd) const
    {
        return axis_sign * bd.adjust_direction(pos - this->origin_)[axis_index];
    }

    //XXX take care. the actual force that would be applied to a particle is
    //    `-dV/dx * calc_force_direction()`.
    coordinate_type calc_force_direction(
            const coordinate_type& pos, const boundary_type& bd) const
    {
        return normal_axis_type::invoke();
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

    real_type cutoff_, margin_, current_margin_;
    coordinate_type origin_; // representative position of a plane.
    std::vector<std::size_t> neighbors_;   // being inside of cutoff range
    std::vector<std::size_t> participant_; // particle that interacts with
};
template<typename traitsT, template<typename> class NormalAxis>
constexpr std::size_t AxisAlignedPlane<traitsT, NormalAxis>::axis_index;
template<typename traitsT, template<typename> class NormalAxis>
constexpr typename AxisAlignedPlane<traitsT, NormalAxis>::real_type
    AxisAlignedPlane<traitsT, NormalAxis>::axis_sign;

template<typename traitsT, template<typename> class NormalAxis>
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

template<typename traitsT, template<typename> class NormalAxis>
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
