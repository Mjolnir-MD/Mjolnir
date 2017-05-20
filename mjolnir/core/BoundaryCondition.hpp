#ifndef MJOLNIR_BOUNDARY_CONDITION
#define MJOLNIR_BOUNDARY_CONDITION
#include <cassert>

namespace mjolnir
{

template<typename realT, typename coordT>
struct UnlimitedBoundary
{
    typedef realT  real_type;
    typedef coordT coordiante_type;

    UnlimitedBoundary() = default;
    ~UnlimitedBoundary() = default;

    coordinate_type adjust_direction(coordinate_type dr) const {return dr;}
    coordinate_type adjust_position (coordinate_type r ) const {return r;}
};

template<typename realT, typename coordT>
struct CubicPeriodicBoundary
{
  public:
    typedef realT  real_type;
    typedef coordT coordiante_type;

  public:
    CubicPeriodicBoundary() = default;
    ~CubicPeriodicBoundary() = default;
    CubicPeriodicBoundary(const coordinate_type& lw, const coordinate_type& up)
    : lower(lw), upper(up), system_size(up-lw), system_size_half(0.5*(up-lw))
    {}

    coordinate_type adjust_direction(coordinate_type dr) const
    {
             if(dr[0] < -system_size_half[0]) dr[0] += system_size[0];
        else if(dr[0] >  system_size_half[0]) dr[0] -= system_size[0];
             if(dr[1] < -system_size_half[1]) dr[1] += system_size[1];
        else if(dr[1] >  system_size_half[1]) dr[1] -= system_size[1];
             if(dr[2] < -system_size_half[2]) dr[2] += system_size[2];
        else if(dr[2] >  system_size_half[2]) dr[2] -= system_size[2];
        return dr;
    }

    coordinate_type adjust_position(coordinate_type pos) const
    {
             if(pos[0] < lower[0]) pos[0] += system_size_[0];
        else if(pos[0] > upper[0]) pos[0] -= system_size_[0];
             if(pos[1] < lower[1]) pos[1] += system_size_[1];
        else if(pos[1] > upper[1]) pos[1] -= system_size_[1];
             if(pos[2] < lower[2]) pos[2] += system_size_[2];
        else if(pos[2] > upper[2]) pos[2] -= system_size_[2];
        return pos;
    }

    coordinate_type&       lower_bound()       {return lower;}
    coordinate_type const& lower_bound() const {return lower;}
    coordinate_type&       upper_bound()       {return upper;}
    coordinate_type const& upper_bound() const {return upper;}
    coordinate_type const& range() const {return system_size;}

  private:

    coordinate_type lower;
    coordinate_type upper;
    coordinate_type system_size;
    coordinate_type system_size_half;
};

}//mjolnir
#endif /* MJOLNIR_BOUNDARY_CONDITION */
