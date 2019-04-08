#ifndef MJOLNIR_CORE_BOUNDARY_CONDITION_HPP
#define MJOLNIR_CORE_BOUNDARY_CONDITION_HPP
#include <mjolnir/math/math.hpp>
#include <cstddef>

namespace mjolnir
{

template<typename realT, typename coordT>
struct UnlimitedBoundary
{
  public:
    using real_type       = realT;
    using coordinate_type = coordT;

  public:
    UnlimitedBoundary() = default;
    ~UnlimitedBoundary() = default;

    coordinate_type adjust_direction(coordinate_type dr) const noexcept {return dr;}
    coordinate_type adjust_position (coordinate_type r ) const noexcept {return r;}
};

template<typename realT, typename coordT>
struct CuboidalPeriodicBoundary
{
  public:
    using real_type       = realT;
    using coordinate_type = coordT;

  public:

    CuboidalPeriodicBoundary(
            const coordinate_type& lw, const coordinate_type& up) noexcept
        : lower_(lw), upper_(up), width_(up - lw), halfw_((up - lw) / 2)
    {}
    ~CuboidalPeriodicBoundary() = default;

    coordinate_type adjust_direction(coordinate_type dr) const noexcept
    {
        if     (math::X(dr) < -math::X(halfw_)) {math::X(dr) += math::X(width_);}
        else if(math::X(dr) >= math::X(halfw_)) {math::X(dr) -= math::X(width_);}
        if     (math::Y(dr) < -math::Y(halfw_)) {math::Y(dr) += math::Y(width_);}
        else if(math::Y(dr) >= math::Y(halfw_)) {math::Y(dr) -= math::Y(width_);}
        if     (math::Z(dr) < -math::Z(halfw_)) {math::Z(dr) += math::Z(width_);}
        else if(math::Z(dr) >= math::Z(halfw_)) {math::Z(dr) -= math::Z(width_);}
        return dr;
    }

    coordinate_type adjust_position(coordinate_type pos) const noexcept
    {
        if     (math::X(pos) <  math::X(lower_)) {math::X(pos) += math::X(width_);}
        else if(math::X(pos) >= math::X(upper_)) {math::X(pos) -= math::X(width_);}
        if     (math::Y(pos) <  math::Y(lower_)) {math::Y(pos) += math::Y(width_);}
        else if(math::Y(pos) >= math::Y(upper_)) {math::Y(pos) -= math::Y(width_);}
        if     (math::Z(pos) <  math::Z(lower_)) {math::Z(pos) += math::Z(width_);}
        else if(math::Z(pos) >= math::Z(upper_)) {math::Z(pos) -= math::Z(width_);}
        return pos;
    }

    coordinate_type const& lower_bound() const noexcept {return lower_;}
    coordinate_type const& upper_bound() const noexcept {return upper_;}
    coordinate_type const& width()       const noexcept {return width_;}

    void set_boundary(const coordinate_type& lw,
                      const coordinate_type& up) noexcept
    {
        this->lower_ = lw;
        this->upper_ = up;
        this->width_ = this->upper_ - this->lower_;
        this->halfw_ = this->width_ * 0.5;
        return;
    }

    void set_lower_bound(const coordinate_type& lw) noexcept
    {
        this->lower_ = lw;
        this->width_ = this->upper_ - this->lower_;
        this->halfw_ = this->width_ * 0.5;
        return;
    }
    void set_upper_bound(const coordinate_type& up) noexcept
    {
        this->upper_ = up;
        this->width_ = this->upper_ - this->lower_;
        this->halfw_ = this->width_ * 0.5;
        return;
    }

  private:

    coordinate_type lower_;
    coordinate_type upper_;
    coordinate_type width_;
    coordinate_type halfw_;
};

}//mjolnir
#endif /* MJOLNIR_BOUNDARY_CONDITION */
