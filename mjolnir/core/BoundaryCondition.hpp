#ifndef MJOLNIR_BOUNDARY_CONDITION
#define MJOLNIR_BOUNDARY_CONDITION
#include <cassert>

namespace mjolnir
{

template<typename traitsT>
struct UnlimitedBoundary
{
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

    UnlimitedBoundary() = default;
    ~UnlimitedBoundary() = default;

    static coordinate_type
    adjust_direction(coordinate_type dr) {return dr;}
    static coordinate_type
    adjust_absolute(coordinate_type r) {return r;}
};

template<typename traitsT>
struct PeriodicBoundaryXYZ
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:
    PeriodicBoundaryXYZ() = default;
    ~PeriodicBoundaryXYZ() = default;

    static void
    set_system(const coordinate_type& lower, const coordinate_type upper)
    {
        lower_            = lower;
        upper_            = upper;
        sys_size_         = upper_ - lower_;
        sys_size_half_pos =  0.5 * sys_size_;
        sys_size_half_neg = -0.5 * sys_size_;
        return;
    }

    static coordinate_type
    adjust_direction(coordinate_type dr)
    {
             if(dr[0] < sys_size_half_neg[0]) dr[0] += sys_size_[0];
        else if(dr[0] > sys_size_half_pos[0]) dr[0] -= sys_size_[0];
             if(dr[1] < sys_size_half_neg[1]) dr[1] += sys_size_[1];
        else if(dr[1] > sys_size_half_pos[1]) dr[1] -= sys_size_[1];
             if(dr[2] < sys_size_half_neg[2]) dr[2] += sys_size_[2];
        else if(dr[2] > sys_size_half_pos[2]) dr[2] -= sys_size_[2];
        return dr;
    }

    static coordinate_type
    adjust_absolute(coordinate_type pos)
    {
             if(pos[0] < lower_[0]) pos[0] += sys_size_[0];
        else if(pos[0] > upper_[0]) pos[0] -= sys_size_[0];
             if(pos[1] < lower_[1]) pos[1] += sys_size_[1];
        else if(pos[1] > upper_[1]) pos[1] -= sys_size_[1];
             if(pos[2] < lower_[2]) pos[2] += sys_size_[2];
        else if(pos[2] > upper_[2]) pos[2] -= sys_size_[2];
        return pos;
    }

  private:

    static coordinate_type lower_;
    static coordinate_type upper_;
    static coordinate_type sys_size_;
    static coordinate_type sys_size_half_pos;
    static coordinate_type sys_size_half_neg;
};

template<typename traitsT>
typename PeriodicBoundaryXYZ<traitsT>::coordinate_type
PeriodicBoundaryXYZ<traitsT>::lower_;
template<typename traitsT>
typename PeriodicBoundaryXYZ<traitsT>::coordinate_type
PeriodicBoundaryXYZ<traitsT>::upper_;
template<typename traitsT>
typename PeriodicBoundaryXYZ<traitsT>::coordinate_type
PeriodicBoundaryXYZ<traitsT>::sys_size_;
template<typename traitsT>
typename PeriodicBoundaryXYZ<traitsT>::coordinate_type
PeriodicBoundaryXYZ<traitsT>::sys_size_half_pos;
template<typename traitsT>
typename PeriodicBoundaryXYZ<traitsT>::coordinate_type
PeriodicBoundaryXYZ<traitsT>::sys_size_half_neg;


}//mjolnir
#endif /* MJOLNIR_BOUNDARY_CONDITION */
