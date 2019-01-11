#ifndef MJOLNIR_BOUNDARY_CONDITION
#define MJOLNIR_BOUNDARY_CONDITION
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
        if     (dr[0] <  -halfw_[0]) {dr[0] += width_[0];}
        else if(dr[0] >=  halfw_[0]) {dr[0] -= width_[0];}
        if     (dr[1] <  -halfw_[1]) {dr[1] += width_[1];}
        else if(dr[1] >=  halfw_[1]) {dr[1] -= width_[1];}
        if     (dr[2] <  -halfw_[2]) {dr[2] += width_[2];}
        else if(dr[2] >=  halfw_[2]) {dr[2] -= width_[2];}
        return dr;
    }

    coordinate_type adjust_position(coordinate_type pos) const noexcept
    {
        if     (pos[0] <  lower_[0]) {pos[0] += width_[0];}
        else if(pos[0] >= upper_[0]) {pos[0] -= width_[0];}
        if     (pos[1] <  lower_[1]) {pos[1] += width_[1];}
        else if(pos[1] >= upper_[1]) {pos[1] -= width_[1];}
        if     (pos[2] <  lower_[2]) {pos[2] += width_[2];}
        else if(pos[2] >= upper_[2]) {pos[2] -= width_[2];}
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
