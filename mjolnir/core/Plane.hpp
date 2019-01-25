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
    using traits_type     = traitsT;
    using system_type     = System<traits_type>;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using boundary_type   = typename traits_type::boundary_type;

    static_assert(!std::is_same<boundary_type,
            CuboidalPeriodicBoundary<real_type, coordinate_type>>::value,
    "Plane that is NOT aligned with axes is not suitable to PeriodicBoundary");

  public:

    Plane(const coordinate_type& pos, const coordinate_type& n,
          const real_type margin = 1)
        : position_(pos), normal_(n), margin_(margin), current_margin_(-1)
    {}

    //XXX   can be negative! because normal vector can define the direction...
    //TODO: consider more clear name
    real_type calc_distance(
            const coordinate_type& pos, const boundary_type& bd) const
    {
        return math::dot_product(this->normal_,
                bd.adjust_direction(pos - this->position_));
    }

    coordinate_type calc_force_direction(
            const coordinate_type& pos, const boundary_type& bd) const
    {
        // calculate `f` that will be used in this form `F = -dV * f`;
        const real_type sign = std::copysign(real_type(1.0), math::dot_product(
                this->normal_, bd.adjust_direction(pos - this->position_)));

        return sign * normal_;
    }

    template<typename Potential>
    void initialize(const system_type& sys, const Potential& pot)
    {
        // update potential parameter
        this->cutoff_      = pot.max_cutoff_length();
        this->participant_ = pot.participants();
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

    void make  (const system_type& sys);
    void update(const real_type dm, const system_type& sys);

    std::vector<std::size_t> const& neighbors() const noexcept
    {return this->negihbors_;}

  private:

    coordinate_type position_; // representative position.
    coordinate_type normal_;   // normal vector

    real_type cutoff_, margin_, current_margin_;
    std::vector<std::size_t> neighbors_;   // being inside of cutoff range
    std::vector<std::size_t> participant_; // particle that interacts with
};

template<typename traitsT>
void Plane<traitsT>::update(const real_type dmargin, const system_type& sys)
{
    this->current_margin_ -= dmargin;
    if(this->current_margin_ < 0)
    {
        this->make(sys);
    }
    return;
}

template<typename traitsT>
void Plane<traitsT>::make(const system_type& sys)
{
    this->neighbors_.clear();
    const real_type threshold = this->cutoff_ * (1 + this->margin_);

    for(std::size_t i : this->participant_)
    {
        const auto d = this->calc_distance(sys[i].position, sys.boundary());
        if(dist < threshold)
        {
            this->neighbors_.push_back(i);
        }
    }

    this->current_margin_ = this->cutoff_ * this->margin_;
    return;
}

} // mjolnir
#endif // MJOLNIR_PLANE_HPP
