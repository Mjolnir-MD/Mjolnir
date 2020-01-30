#ifndef MJOLNIR_FORCEFIELD_AFMFIT_AFMFIT_INTEARACTION_HPP
#define MJOLNIR_FORCEFIELD_AFMFIT_AFMFIT_INTEARACTION_HPP
#include <mjolnir/math/functions.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/core/ExternalForceInteractionBase.hpp>
#include <utility>
#include <algorithm>
#include <vector>

namespace mjolnir
{

// AFM Fitting Interaction.
//
// This is an implementation of AFM Flexible Fitting Potential developed by
// the following paper.
// - Toru Niina, Sotaro Fuchigami, Shoji Takada (2020) JCTC
// - doi:10.1021/acs.jctc.9b00991
//
// The functional form of the potential is too complicated to write here.
//
// Since this interaction tightly coupled with the potential function, it is
// implemented as a single class.
//
// XXX: This implementation is based on the following assumptions.
// - Image is recutangular shape
// - Each pixel shows the height at the center of pixel region.
// - The origin of the image is (0, 0). That means that the first pixel shows
//   the height at (pixel_size_x/2, pixel_size_y/2).
// - The image contains pixels in the following order.
//   - (x0, y0), (x1, y0), (x2, y0), ... (xn, y0),
//     (x0, y1), (x1, y1), (x2, y1), ... (xn, y1),
//     ...
//     (x0, ym), (x1, ym), (x2, ym), ... (xn, ym),
//
template<typename traitsT>
class AFMFitInteraction final : public ExternalForceInteractionBase<traitsT>
{
  public:
    using traits_type     = traitsT;
    using base_type       = ExternalForceInteractionBase<traits_type>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using boundary_type   = typename base_type::boundary_type;
    using parameter_type  = real_type; // radius of the particle

    // radius parameter
    static constexpr real_type default_parameter() {return 0.0;}

  public:

    // It consider the first pixel locates (pixel_size_x/2, pixel_size_y/2).
    // Take care when you visualize it.
    AFMFitInteraction(const real_type k, const real_type gamma,
        const real_type z0, const real_type cutoff, const real_type margin,
        const real_type   sgm_x, const real_type   sgm_y,
        const real_type   pix_x, const real_type   pix_y,
        const std::size_t len_x, const std::size_t len_y,
        const std::vector<std::pair<std::size_t, parameter_type>>& params,
        const std::vector<real_type>& image)
        : k_(k), gamma_(gamma), rgamma_(1 / gamma), z0_(z0),
          stage_term_(std::exp(z0 / gamma)), cutoff_(cutoff), margin_(margin),
          sgm_x_ (sgm_x), sgm_y_ (sgm_y), rsgm_x_(1 / sgm_x), rsgm_y_(1 / sgm_y),
          rsgm_x_sq_(rsgm_x_ * rsgm_x_),  rsgm_y_sq_(rsgm_y_ * rsgm_y_),
          pxl_x_ (pix_x), pxl_y_ (pix_y), rpxl_x_(1 / pix_x), rpxl_y_(1 / pix_y),
          len_x_ (len_x), len_y_ (len_y),
          pixel_in_cutoff_x_(std::ceil(sgm_x_ * cutoff_ * (1 + margin_) * rpxl_x_)),
          pixel_in_cutoff_y_(std::ceil(sgm_y_ * cutoff_ * (1 + margin_) * rpxl_y_)),
          max_pixel_interact_((1 + 2 * pixel_in_cutoff_x_) *
                              (1 + 2 * pixel_in_cutoff_y_)),
          rH_ref_ref_sum_(real_type(1) / std::accumulate(
              image.begin(), image.end(), real_type(0),
              [](const real_type init, const real_type H) noexcept -> real_type {
                  return init + H * H;
              })),
          H_ref_ (image), pixel_positions_(image.size()), H_sim_(image.size(), 0),
          sumexp_(image.size(), 0), rsumexp_(image.size(), 0),
          rH_ref_sim_sum_(0), rH_sim_sim_sum_(0), current_margin_(-1),
          particle_to_pixel_(max_pixel_interact_ * params.size(),
                             std::numeric_limits<std::size_t>::max())
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("z0 = ", this->z0_, "stage term = ", this->stage_term_);
        MJOLNIR_LOG_INFO(1+2*pixel_in_cutoff_x_, 'x', 1+2*pixel_in_cutoff_y_,
                " pixels exist in the cutoff range (", cutoff_, " sigma)");

        // initialize participants and their radius parameters
        this->parameters_  .resize (params.size(), default_parameter());
        this->participants_.reserve(params.size());
        for(const auto& idxp : params)
        {
            const auto idx = idxp.first;
            this->participants_.push_back(idx);
            if(idx >= this->parameters_.size())
            {
                this->parameters_.resize(idx+1, default_parameter());
            }
            this->parameters_.at(idx) = idxp.second;
        }

        // intiialize pixel positions
        for(std::size_t yi=0; yi < len_y_; ++yi)
        {
            const std::size_t offset =  yi        * len_x_;
            const real_type   y0     = (yi + 0.5) * pxl_y_;

            for(std::size_t xi=0; xi < len_x_; ++xi)
            {
                const real_type x0 = (xi + 0.5) * pxl_x_;
                this->pixel_positions_.at(offset + xi) = std::make_pair(x0, y0);
                MJOLNIR_LOG_INFO("position of (", xi, ", ", yi, ") th pixel is "
                                 "(", x0, ", ", y0, ")");
            }
        }
    }
    ~AFMFitInteraction() override {}

    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());

        this->update(sys);
        return;
    }

    void update(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());
        this->construct_list(sys);
        return;
    }

    void reduce_margin(const real_type dmargin, const system_type& sys) override
    {
        this->current_margin_ -= dmargin;
        if(this->current_margin_ < 0)
        {
            this->construct_list(sys);
        }
        return;
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        const auto abs_cutoff_x = this->cutoff_ * sgm_x_;
        const auto abs_cutoff_y = this->cutoff_ * sgm_y_;
        const auto scaled_margin_x = (abs_cutoff_x + current_margin_) * scale - abs_cutoff_x;
        const auto scaled_margin_y = (abs_cutoff_y + current_margin_) * scale - abs_cutoff_y;

        this->current_margin_ = std::min(scaled_margin_x, scaled_margin_y);
        if(this->current_margin_ < 0)
        {
            this->construct_list(sys);
        }
        return;
    }

    void      calc_force (system_type&)       const noexcept override;
    real_type calc_energy(system_type const&) const noexcept override;

    std::string name() const override {return "AFMFlexibleFitting";}

    base_type* clone() const override
    {
        std::vector<std::pair<std::size_t, parameter_type>> params;
        params.reserve(participants_.size());
        for(const auto i : this->participants_)
        {
            params.emplace_back(i, this->parameters_.at(i));
        }
        return new AFMFitInteraction(k_, gamma_, z0_, cutoff_, margin_,
            sgm_x_, sgm_y_, pxl_x_, pxl_y_, len_x_, len_y_, params, H_ref_);
    }

    // accessor for testing.

    real_type   k()        const noexcept {return k_;}
    real_type   gamma()    const noexcept {return gamma_;}
    real_type   z0()       const noexcept {return z0_;}
    real_type   cutoff()   const noexcept {return cutoff_;}
    real_type   margin()   const noexcept {return margin_;}
    real_type   sigma_x()  const noexcept {return sgm_x_;}
    real_type   sigma_y()  const noexcept {return sgm_y_;}
    real_type   pixel_x()  const noexcept {return pxl_x_;}
    real_type   pixel_y()  const noexcept {return pxl_y_;}
    std::size_t length_x() const noexcept {return len_x_;}
    std::size_t length_y() const noexcept {return len_y_;}

    std::vector<std::size_t> const& participants() const noexcept {return participants_;}
    std::vector<real_type>   const& parameters()   const noexcept {return parameters_;}
    std::vector<real_type>   const& image()        const noexcept {return H_ref_;}

  private:

    // it calculates correlation coefficient and update pixel list
    real_type calc_correlation(const system_type& sys) const;

    void construct_list(const system_type& sys);

  private:

    real_type   k_, gamma_, rgamma_, z0_, stage_term_, cutoff_, margin_;
    real_type   sgm_x_,  sgm_y_, rsgm_x_, rsgm_y_, rsgm_x_sq_, rsgm_y_sq_;
    real_type   pxl_x_,  pxl_y_, rpxl_x_, rpxl_y_;
    std::size_t len_x_, len_y_;
    std::size_t pixel_in_cutoff_x_, pixel_in_cutoff_y_, max_pixel_interact_;

    // particle radius
    std::vector<parameter_type> parameters_;
    std::vector<std::size_t>    participants_;

    // reference image
    real_type   rH_ref_ref_sum_;
    std::vector<real_type> H_ref_;
    std::vector<std::pair<real_type, real_type>> pixel_positions_;

    // generated image and the intermediate results.
    // These are marked mutable because calc_correlation modifies these.
    // These are large, dynamically-allocated vector. Allocating this each
    // `calc_force` decreases the runtime efficiency.
    mutable std::vector<real_type> H_sim_;
    mutable std::vector<real_type> sumexp_, rsumexp_;
    mutable real_type rH_ref_sim_sum_, rH_sim_sim_sum_;

    // A kind of cell list. It manages which particle corresponds to which pixel.
    // pid -> {pxid, pxid, ...}
    // {px1_pid1, px2_pid1, ..., pxN_pid1, px1_pid2, ... pxN_pidM}
    real_type current_margin_;
    std::vector<std::size_t> particle_to_pixel_;
};

template<typename traitsT>
void AFMFitInteraction<traitsT>::construct_list(const system_type& sys)
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();
    assert(this->particle_to_pixel_.size() ==
           this->participants_.size() * this->max_pixel_interact_);

    std::fill(particle_to_pixel_.begin(), particle_to_pixel_.end(),
              std::numeric_limits<std::size_t>::max());

    const auto x_cutoff_margin = sgm_x_ * cutoff_ * (1 + margin_);
    const auto y_cutoff_margin = sgm_y_ * cutoff_ * (1 + margin_);

    // the edge of the image
    const int max_x_img = static_cast<int>(len_x_ - 1);
    const int max_y_img = static_cast<int>(len_y_ - 1);
    std::size_t participants_idx = 0;
    for(const std::size_t i : this->participants_)
    {
        const auto& pos = sys.position(i);

        const auto pos_x = math::X(pos);
        const auto pos_y = math::Y(pos);

        // center position in the pixel number coordinate
        const int ctr_xi = static_cast<int>(std::ceil(pos_x * rpxl_x_));
        const int ctr_yi = static_cast<int>(std::ceil(pos_y * rpxl_y_));

        // related pixel region. clamped in the region, [0, max_index]
        const int min_xi = math::clamp<int>(ctr_xi - pixel_in_cutoff_x_, 0, max_x_img);
        const int max_xi = math::clamp<int>(ctr_xi + pixel_in_cutoff_x_, 0, max_x_img);
        const int min_yi = math::clamp<int>(ctr_yi - pixel_in_cutoff_y_, 0, max_y_img);
        const int max_yi = math::clamp<int>(ctr_yi + pixel_in_cutoff_y_, 0, max_y_img);

        std::size_t pxlist_offset = participants_idx * this->max_pixel_interact_;

        for(int yi = min_yi; yi <= max_yi; ++yi)
        {
            for(int xi = min_xi; xi <= max_xi; ++xi)
            {
                const std::size_t pxl = yi * len_x_ + xi;

                const real_type x0 = this->pixel_positions_[pxl].first;
                const real_type y0 = this->pixel_positions_[pxl].second;

                const real_type dx = pos_x - x0;
                const real_type dy = pos_y - y0;

                if(std::abs(dx) <= x_cutoff_margin &&
                   std::abs(dy) <= y_cutoff_margin)
                {
                    particle_to_pixel_[pxlist_offset] = pxl;
                    pxlist_offset += 1;
                }
            }
        }
        MJOLNIR_LOG_DEBUG("position (", pos_x, ", ", pos_y, ") relates ",
            pxlist_offset - participants_idx * this->max_pixel_interact_,
            " pixels.");

        assert(pxlist_offset <= (participants_idx + 1) * this->max_pixel_interact_);
        participants_idx += 1;
    }
    assert(participants_idx == participants_.size());
    this->current_margin_ = std::min(sgm_x_ * cutoff_ * margin_,
                                     sgm_y_ * cutoff_ * margin_);
    return;
}

template<typename traitsT>
typename AFMFitInteraction<traitsT>::real_type
AFMFitInteraction<traitsT>::calc_correlation(const system_type& sys) const
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();
    assert(this->H_sim_.size()           == this->H_ref_.size());
    assert(this->rsumexp_.size()         == this->H_ref_.size());
    assert(this->pixel_positions_.size() == this->H_ref_.size());

    std::fill(this->sumexp_.begin(), this->sumexp_.end(), this->stage_term_);

    const real_type cutoff_x = sgm_x_ * cutoff_;
    const real_type cutoff_y = sgm_y_ * cutoff_;

    std::size_t participants_idx = 0;
    for(const std::size_t i : this->participants_)
    {
        const auto& pos = sys.position(i);
        const auto  rad = parameters_[i];

        const std::size_t pxlist_offset = participants_idx * max_pixel_interact_;
        for(std::size_t px=0; px<this->max_pixel_interact_; ++px)
        {
            const auto pxl = particle_to_pixel_[pxlist_offset + px];
            if(pxl == std::numeric_limits<std::size_t>::max())
            {
                break;
            }

            const real_type x0 = this->pixel_positions_[pxl].first;
            const real_type y0 = this->pixel_positions_[pxl].second;

            const real_type dx = math::X(pos) - x0;
            const real_type dy = math::Y(pos) - y0;
            const real_type dz = math::Z(pos) + rad - z0_;

            if(cutoff_x < std::abs(dx) || cutoff_y < std::abs(dy))
            {
                continue;
            }

            const real_type gx = -dx * dx * rsgm_x_sq_ * real_type(0.5);
            const real_type gy = -dy * dy * rsgm_y_sq_ * real_type(0.5);
            const real_type hz = dz * this->rgamma_;

            sumexp_[pxl] += std::exp(gx + gy + hz);
        }
        participants_idx += 1;
    }
    assert(participants_idx == participants_.size());

    // sumexp_ calculated.

    assert(H_sim_.size() == H_ref_.size());

    real_type H_sim_sim_sum = 0;
    real_type H_ref_sim_sum = 0;
    for(std::size_t pxl = 0; pxl < this->H_sim_.size(); ++pxl)
    {
        MJOLNIR_LOG_DEBUG("sum of exp(x, y, z) = ", sumexp_[pxl]);
        const real_type H = this->z0_ + this->gamma_ * std::log(sumexp_[pxl]);

        this->rsumexp_[pxl] = real_type(1) / this->sumexp_[pxl];
        this->H_sim_[pxl]   = H;

        H_sim_sim_sum += H * H;
        H_ref_sim_sum += this->H_ref_[pxl] * H;

        MJOLNIR_LOG_DEBUG("H_ref = ", H_ref_.at(pxl));
        MJOLNIR_LOG_DEBUG("H_sim = ", H);
    }

    this->rH_sim_sim_sum_ = real_type(1) / H_sim_sim_sum;
    this->rH_ref_sim_sum_ = real_type(1) / H_ref_sim_sum;

    MJOLNIR_LOG_DEBUG("sum of Hsim * Hsim = ", H_sim_sim_sum);
    MJOLNIR_LOG_DEBUG("sum of Href * Hsim = ", H_ref_sim_sum);
    MJOLNIR_LOG_DEBUG("sum of Href * Href = ", 1.0 / rH_ref_ref_sum_);

    return H_ref_sim_sum * std::sqrt(rH_ref_ref_sum_ * rH_sim_sim_sum_);
}

template<typename traitsT>
void AFMFitInteraction<traitsT>::calc_force(system_type& sys) const noexcept
{
    const real_type cc    = this->calc_correlation(sys);
    const real_type coeff = this->k_ * cc;

    const real_type cutoff_x = sgm_x_ * cutoff_;
    const real_type cutoff_y = sgm_y_ * cutoff_;

    std::size_t participants_idx = 0;
    for(const std::size_t i : this->participants_)
    {
        const auto& pos = sys.position(i);
        const auto  rad = parameters_[i];
        auto f = math::make_coordinate<coordinate_type>(0, 0, 0);

        const std::size_t offset = this->max_pixel_interact_ * participants_idx;
        for(std::size_t j=0; j<this->max_pixel_interact_; ++j)
        {
            const std::size_t pxl = this->particle_to_pixel_[offset + j];

            if(pxl == std::numeric_limits<std::size_t>::max())
            {
                break;
            }
            const real_type x0 = this->pixel_positions_[pxl].first;
            const real_type y0 = this->pixel_positions_[pxl].second;

            const real_type dx = math::X(pos) - x0;
            const real_type dy = math::Y(pos) - y0;
            const real_type dz = math::Z(pos) + rad - z0_;

            if(cutoff_x < std::abs(dx) || cutoff_y < std::abs(dy))
            {
                continue;
            }

            const real_type gx = real_type(-0.5) * dx * dx * rsgm_x_sq_;
            const real_type gy = real_type(-0.5) * dy * dy * rsgm_y_sq_;
            const real_type hz = dz * this->rgamma_;

            const real_type term = coeff * std::exp(gx + gy + hz) *
                (this->H_ref_[pxl] * this->rH_ref_sim_sum_ -
                 this->H_sim_[pxl] * this->rH_sim_sim_sum_) *
                 this->gamma_ * this->rsumexp_[pxl];

            math::X(f) += term * (-dx * this->rsgm_x_sq_);
            math::Y(f) += term * (-dy * this->rsgm_y_sq_);
            math::Z(f) += term * this->rgamma_;
        }

        sys.force(i) += f;
        participants_idx += 1;
    }
    assert(participants_idx == participants_.size());

    sys.attribute("correlation") = cc;
    return;
}

template<typename traitsT>
typename AFMFitInteraction<traitsT>::real_type
AFMFitInteraction<traitsT>::calc_energy(const system_type& sys) const noexcept
{
    return this->k_ * (real_type(1) - this->calc_correlation(sys));
}

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize major use-cases
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

namespace mjolnir
{

extern template class AFMFitInteraction<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class AFMFitInteraction<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class AFMFitInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class AFMFitInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif//MJOLNIR_AFM_INTEARACTION_BASE
