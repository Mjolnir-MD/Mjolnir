#ifndef MJOLNIR_CORE_AXIS_ALIGNED_PLANE_HPP
#define MJOLNIR_CORE_AXIS_ALIGNED_PLANE_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

/* These classes represents Normal Vectors. (+/- * XYZ) */

template<typename realT>
class PositiveDirection
{
  public:
    using real_type = realT;
    static constexpr real_type sign() {return real_type(1);}
};
template<typename realT>
class NegativeDirection
{
  public:
    using real_type = realT;
    static constexpr real_type sign() {return real_type(-1);}
};

template<typename traitsT, template<typename> class Direction>
class XAxis
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using direction_type  = Direction<real_type>;

    static constexpr real_type sign() noexcept {return direction_type::sign();}

    static coordinate_type get_normal_vector() noexcept
    {
        return math::make_coordinate<coordinate_type>(
                direction_type::sign() * 1.0, 0, 0);
    }
    static real_type& get_coordinate(coordinate_type&       v) noexcept {return math::X(v);}
    static real_type  get_coordinate(coordinate_type const& v) noexcept {return math::X(v);}
};

template<typename traitsT, template<typename> class Direction>
class YAxis
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using direction_type  = Direction<real_type>;

    static constexpr real_type sign() noexcept {return direction_type::sign();}

    static coordinate_type get_normal_vector() noexcept
    {
        return math::make_coordinate<coordinate_type>(
                0, direction_type::sign() * 1.0, 0);
    }
    static real_type& get_coordinate(coordinate_type&       v) noexcept {return math::Y(v);}
    static real_type  get_coordinate(coordinate_type const& v) noexcept {return math::Y(v);}
};

template<typename traitsT, template<typename> class Direction>
class ZAxis
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using direction_type  = Direction<real_type>;

    static constexpr real_type sign() noexcept {return direction_type::sign();}

    static coordinate_type get_normal_vector() noexcept
    {
        return math::make_coordinate<coordinate_type>(
                0, 0, direction_type::sign() * 1.0);
    }
    static real_type& get_coordinate(coordinate_type&       v) noexcept {return math::Z(v);}
    static real_type  get_coordinate(coordinate_type const& v) noexcept {return math::Z(v);}
};

template<typename traitsT>
using PositiveXDirection = XAxis<traitsT, PositiveDirection>;
template<typename traitsT>
using NegativeXDirection = XAxis<traitsT, NegativeDirection>;

template<typename traitsT>
using PositiveYDirection = YAxis<traitsT, PositiveDirection>;
template<typename traitsT>
using NegativeYDirection = YAxis<traitsT, NegativeDirection>;

template<typename traitsT>
using PositiveZDirection = ZAxis<traitsT, PositiveDirection>;
template<typename traitsT>
using NegativeZDirection = ZAxis<traitsT, NegativeDirection>;

/*! @brief axis aligned Plane for ExternalDistanceInteraction.                *
 *  @details represents flat Plane. It provides a method to calculate         *
 *           a distance between particle and the plane, force direction       *
 *           a position of a particle. It also provide a functionality of     *
 *           a neighbor-list.                                                 */
template<typename traitsT, typename axisT>
class AxisAlignedPlane
{
  public:
    using traits_type      = traitsT;
    using axis_type        = axisT;
    using system_type      = System<traits_type>;
    using real_type        = typename traits_type::real_type;
    using coordinate_type  = typename traits_type::coordinate_type;
    using boundary_type    = typename traits_type::boundary_type;

  public:

    AxisAlignedPlane(const real_type position, const real_type margin = 1)
        : margin_(margin), current_margin_(-1),
          origin_(math::make_coordinate<coordinate_type>(0,0,0))
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        axis_type::get_coordinate(this->origin_) = position;

        MJOLNIR_LOG_INFO("position        = ", position);
        MJOLNIR_LOG_INFO("origin of plane = ", this->origin_);
        MJOLNIR_LOG_INFO("margin          = ", this->margin_);
    }

    //XXX can be negative! because normal vector can define the direction...
    inline real_type
    calc_distance(const coordinate_type& pos, const boundary_type& bd) const
    {
        return axis_type::sign() * axis_type::get_coordinate(
                bd.adjust_direction(pos - this->origin_));
    }

    //XXX take care. the actual force that would be applied to a particle is
    //    `-dV/dx * calc_force_direction()`.
    coordinate_type calc_force_direction(
            const coordinate_type&, const boundary_type&) const
    {
        return axis_type::get_normal_vector();
    }

    template<typename Potential>
    void initialize(const system_type& sys, const Potential& pot)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // assuming pot is already initialized!
        this->cutoff_         = pot.max_cutoff_length();
        this->participant_    = pot.participants();

        MJOLNIR_LOG_INFO("cutoff         = ", this->cutoff_);
        MJOLNIR_LOG_INFO("margin         = ", this->margin_);
        MJOLNIR_LOG_INFO("participant_   = ", this->origin_);

        this->make(sys);
        MJOLNIR_LOG_INFO("current margin = ", this->current_margin_);

        return;
    }

    // neighbor-list stuff
    void make  (const system_type& sys)
    {
        this->neighbors_.clear();
        const real_type threshold = this->cutoff_ * (1 + this->margin_);
        const real_type thr2 = threshold * threshold;

        for(std::size_t i : this->participant_)
        {
            const auto d = this->calc_distance(sys.position(i), sys.boundary());
            if(d * d <= thr2)
            {
                this->neighbors_.push_back(i);
            }
        }

        this->current_margin_ = this->cutoff_ * this->margin_;
        return;
    }

    void update(const real_type dmargin, const system_type& sys)
    {
        this->current_margin_ -= dmargin;
        if(this->current_margin_ < 0)
        {
            this->make(sys);
        }
        return;
    }

    std::vector<std::size_t> const& neighbors() const noexcept
    {return this->neighbors_;}

  private:

    real_type cutoff_, margin_, current_margin_;
    coordinate_type origin_; // representative position of a plane.
    std::vector<std::size_t> neighbors_;   // being inside of cutoff range
    std::vector<std::size_t> participant_; // particle that interacts with
};

} // mjolnir
#endif // MJOLNIR_AXIS_ALIGNED_PLANE_HPP
