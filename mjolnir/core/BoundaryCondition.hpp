#ifndef MJOLNIR_BOUNDARY_CONDITION
#define MJOLNIR_BOUNDARY_CONDITION

namespace mjolnir
{

template<typename traitsT>
struct UnlimitedBoundary
{
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

    UnlimitedBoundary() = default;
    UnlimitedBoundary(const coordinate_type&){}
    ~UnlimitedBoundary() = default;

    coordinate_type
    operator()(coordinate_type dr) const {return dr;}
};

template<typename traitsT>
struct PeriodicBoundaryXYZ
{
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

    PeriodicBoundaryXYZ() = default;
    ~PeriodicBoundaryXYZ() = default;

    PeriodicBoundaryXYZ(const coordinate_type& sys_size)
        : sys_size_(sys_size),
          sys_size_half_pos(0.5 * sys_size), sys_size_half_neg(-0.5 * sys_size)
    {}

    void set_system_size(const coordinate_type& sys_size)
    {
        this->sys_size_         = sys_size;
        this->sys_size_half_pos = 0.5 * sys_size;
        this->sys_size_half_neg = -0.5 * sys_size;
        return;
    }

    coordinate_type
    operator()(coordinate_type dr) const
    {
             if(dr[0] < sys_size_half_neg[0]) dr[0] += sys_size_[0];
        else if(dr[0] > sys_size_half_pos[0]) dr[0] -= sys_size_[0];
             if(dr[1] < sys_size_half_neg[1]) dr[1] += sys_size_[1];
        else if(dr[1] > sys_size_half_pos[1]) dr[1] -= sys_size_[1];
             if(dr[2] < sys_size_half_neg[2]) dr[2] += sys_size_[2];
        else if(dr[2] > sys_size_half_pos[2]) dr[2] -= sys_size_[2];
        return dr;
    }

  private:
    coordinate_type sys_size_;
    coordinate_type sys_size_half_pos;
    coordinate_type sys_size_half_neg;
};



}//mjolnir
#endif /* MJOLNIR_BOUNDARY_CONDITION */
