#ifndef MJOLNIR_PLANE_HPP
#define MJOLNIR_PLANE_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{

/*! @brief Plane for ExternalDistanceInteraction.                             *
 *  @details represents flat Plane. It provides a method to calculate         *
 *           a distance between particle and the plane, force direction       *
 *           a position of a particle. It also provide a functionality of     *
 *           a neighbor-list.                                                 */
template<typename traitsT>
class Plane
{
  public:
    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type        real_type;
    typedef typename traits_type::coordinate_type  coordinate_type;
    typedef typename traits_type::boundary_type    boundary_type;
    typedef typename traits_type::particle_type    particle_type;

    static_assert(!std::is_same<boundary_type,
            CubicPeriodicBoundary<real_type, coordinate_type>>::value,
    "Plane that is NOT aligned with axes is not suitable to PeriodicBoundary");

  public:

    Plane(const coordinate_type& pos, const coordinate_type& n,
          const real_type mergin = 1)
        : position_(pos), normal_(n), mergin_(mergin), current_mergin_(-1)
    {}

    real_type calc_distance(
            const coordinate_type& pos, const boundary_type& bd) const
    {
        return dot_product(this->normal_,
                bd.adjust_direction(pos - this->position_));
    }

    coordinate_type calc_force_direction(
            const coordinate_type& pos, const boundary_type& bd) const
    {
        // calculate `f` that will be used in this form `F = -dV * f`;
        const real_type sign = std::copysign(real_type(1.0), dot_product(
                this->normal_, bd.adjust_direction(pos - this->position_)));
        return sign * normal_;
    }

    //! after calling it, update() should be called!
    void set_cutoff(const real_type rc) noexcept
    {this->cutoff_ = rc; this->current_mergin_ = -1;}
    void set_mergin(const real_type mg) noexcept
    {this->mergin_ = mg; this->current_mergin_ = -1;}

    template<typename Potential>
    void initialize(const system_type& sys, const Potential& pot)
    {
        const auto& topol = sys.topology();

        // update group ids
        this->group_id_.clear();
        this->group_id_.resize(sys.size());
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            this->group_id_[i] = topol.group_of(i);
        }

        // update potential parameter
        this->ignored_groups_ = pot.ignored_group();    // these functions are
        this->exclusion_      = pot.ignored_particle(); // required by External
        this->cutoff_         = pot.cutoff_length();    // Potential Concept.

        // reconstruct neighbor list
        this->make(sys);
        return;
    }

    void make(const system_type& sys)
    {
        this->neighbors_.clear();
        const real_type threshold = cutoff_ * (1 + mergin_);
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            if(std::find(ignored_groups_.cbegin(), ignored_groups_.cend(),
                         group_id_[i]) != ignored_groups_.cend()) {continue;}
            if(std::find(exclusion_.cbegin(), exclusion_.cend(), i) !=
                         exclusion_.cend()) {continue;}

            const real_type dist =
                this->calc_distance(sys[i].position, sys.boundary());
            if(dist < threshold)
            {
                this->neighbors_.push_back(i);
            }
        }
        this->current_mergin_ = cutoff_ * mergin_;
        return;
    }
    void update(const system_type& sys)
    {
        this->current_mergin_ -= sys.largest_displacement();
        if(this->current_mergin_ < 0)
        {
            this->make(sys);
        }
        return;
    }

    //! update exclusion list and cutoff length.
    template<typename Potential>
    void reconstruct(const system_type& sys, const Potential& pot)
    {
        // do the same thing as the initialization
        this->initialize(sys, pot);
        return;
    }

    std::vector<std::size_t> const& neighbors() const noexcept
    {return this->negihbors_;}

  private:

    coordinate_type position_; // representative position.
    coordinate_type normal_;   // normal vector

    real_type cutoff_, mergin_, current_mergin_;
    std::vector<std::size_t> neighbors_; // particle that interacts with
    std::vector<std::size_t> exclusion_; // particle that should be ignored
    std::vector<std::size_t> group_id_;  // group id of particle (same as topol)
    std::vector<std::size_t> ignored_groups_; // same as pot.ignored_group()
};

} // mjolnir
#endif // MJOLNIR_PLANE_HPP
