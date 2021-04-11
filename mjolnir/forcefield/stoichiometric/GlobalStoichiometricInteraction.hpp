#ifndef MJOLNIR_FORCEFIELD_STOICHIOMETRIC_GLOBAL_STOICHIOMETRIC_INTERACTION_HPP
#define MJOLNIR_FORCEFIELD_STOICHIOMETRIC_GLOBAL_STOICHIOMETRIC_INTERACTION_HPP
#include <mjolnir/forcefield/stoichiometric/GlobalStoichiometricInteractionPotential.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/GlobalInteractionBase.hpp>
#include <mjolnir/core/SpatialPartitionBase.hpp>
#include <mjolnir/util/logger.hpp>
#include <memory>

namespace mjolnir
{

// The interaction which conserve stoichiometry.

template<typename traitsT>
class GlobalStoichiometricInteraction final : public GlobalInteractionBase<traitsT>
{
  public:

    using traits_type            = traitsT;
    using base_type              = GlobalInteractionBase<traitsT>;
    using index_type             = std::size_t;
    using real_type              = typename base_type::real_type;
    using coordinate_type        = typename base_type::coordinate_type;
    using system_type            = typename base_type::system_type;
    using topology_type          = typename base_type::topology_type;
    using potential_type         = GlobalStoichiometricInteractionPotential<traitsT>;
    using partition_type         = SpatialPartition<traitsT, potential_type>;
    using potential_buffer_type  = std::vector<real_type>;
    using derivative_buffer_type = std::vector<coordinate_type>;

    using potential_buffer_iterator    = typename potential_buffer_type::iterator;
    using derivative_buffer_iterator   = typename derivative_buffer_type::iterator;
    using partner_range_type           = typename partition_type::range_type;
    using potential_buffer_range_type  = range<potential_buffer_iterator>;
    using derivative_buffer_range_type = range<derivative_buffer_iterator>;
    using partner_buffer_tuple_type    = std::tuple<partner_range_type,
                                                    potential_buffer_range_type,
                                                    derivative_buffer_range_type>;
    using partner_buffer_ranges_type   = std::vector<partner_buffer_tuple_type>;

  public:
    GlobalStoichiometricInteraction(potential_type&& pot, partition_type&& part,
        real_type epsilon,
        std::size_t coefa, std::size_t coefb)
        : potential_(std::move(pot)), partition_(std::move(part)),
          epsilon_(epsilon), coefa_(coefa), coefb_(coefb)
    {
        const std::size_t participants_a_num = potential_.participants_a_num();
        const std::size_t participants_b_num = potential_.participants_b_num();
        const std::size_t participants_ab    = participants_a_num * participants_b_num;
        derivs_buff_.resize(participants_ab);
        potentials_buff_.resize(participants_ab);
        partner_buffer_ranges_.resize(participants_a_num);

        pots_sum_a_  .resize(participants_a_num);
        pots_sum_b_  .resize(participants_b_num);
        derivs_sum_a_.resize(participants_a_num);
        derivs_sum_b_.resize(participants_b_num);
    }
    ~GlobalStoichiometricInteraction() override {}

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());
        potential_.initialize(sys, topol);
        partition_.initialize(sys, this->potential_);
        this->update_buffer_range();

        // When we calculate potential or derivative for specific idx, we buffering the
        // value in intermediate value container. So, we have to know which index in
        // the system correspond to which index of buffering container. These for loop
        // is to make that container.
        idx_buffer_map_.resize(sys.size());
        const auto participants_a = potential_.participants_a();
        const auto participants_b = potential_.participants_b();
        for(std::size_t idx_a=0; idx_a<participants_a.size(); ++idx_a)
        {
            idx_buffer_map_[participants_a[idx_a]] = idx_a;
        }
        for(std::size_t idx_b=0; idx_b<participants_b.size(); ++idx_b)
        {
            idx_buffer_map_[participants_b[idx_b]] = idx_b;
        }
    }

    void update(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential is ", this->name());
        this->potential_.update(sys,topol);
        this->partition_.initialize(sys, potential_);
        this->update_buffer_range();
    }

    void reduce_margin(const real_type dmargin, const system_type& sys) override
    {
        if(partition_.reduce_margin(dmargin, sys, potential_))
        {
            this->update_buffer_range();
        }
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        partition_.scale_margin(scale, sys, potential_);
    }

    void      calc_force (system_type&)           const noexcept override;
    real_type calc_energy(const system_type&)     const noexcept override;
    real_type calc_force_and_energy(system_type&) const noexcept override;

    std::string name() const override {return "Stoichiometric";}

    potential_type const& potential() const noexcept {return potential_;}
    potential_type &      potential()       noexcept {return potential_;}

    base_type* clone() const override
    {
        return new GlobalStoichiometricInteraction(
                potential_type(potential_), partition_type(partition_),
                epsilon_, coefa_, coefb_);
    }

  private:
    void update_buffer_range()
    {
        MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
        const auto        participants_a     = potential_.participants_a();
        const std::size_t participants_a_num = potential_.participants_a_num();
        auto pot_buff_iter   = potentials_buff_.begin();
        auto deriv_buff_iter = derivs_buff_.begin();

        for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
        {
            const index_type i = participants_a[idx_a];

            const partner_range_type&  partner          = partition_.partners(i);
            const index_type           range_size       = partner.size();
            MJOLNIR_LOG_DEBUG("partner size of particle ", i, " is ", range_size);

            potential_buffer_iterator  pot_first_iter   = pot_buff_iter;
            derivative_buffer_iterator deriv_first_iter = deriv_buff_iter;

            std::advance(pot_buff_iter,   range_size);
            std::advance(deriv_buff_iter, range_size);
            auto pot_range   = make_range(pot_first_iter,   pot_buff_iter);
            auto deriv_range = make_range(deriv_first_iter, deriv_buff_iter);

            partner_buffer_tuple_type& partner_buffer_range = partner_buffer_ranges_[idx_a];
            std::get<0>(partner_buffer_range) = partner;
            std::get<1>(partner_buffer_range) = make_range(pot_first_iter,   pot_buff_iter);
            std::get<2>(partner_buffer_range) = make_range(deriv_first_iter, deriv_buff_iter);
        }
    }

  private:

    potential_type potential_;
    partition_type partition_;

    real_type      epsilon_;
    std::size_t    coefa_;
    std::size_t    coefb_;

    // -----------------------------------------------------------------------
    // Variable for mapping system index to buffering index
    std::vector<std::size_t> idx_buffer_map_;

    // Variables for buffering intermediate value.
    // These value can change in calc_force and calc_energy function.
    mutable potential_buffer_type  potentials_buff_;
    mutable derivative_buffer_type derivs_buff_;
    mutable potential_buffer_type  pots_sum_a_;
    mutable potential_buffer_type  pots_sum_b_;

    // Variables to specify the range of specific a partners in buffer.
    mutable partner_buffer_ranges_type partner_buffer_ranges_;

    // sum of derivation of potential function for specific first particle.
    mutable derivative_buffer_type     derivs_sum_a_;

    // sum of derivation of potential function for specific second particle.
    mutable derivative_buffer_type     derivs_sum_b_;
};

template<typename traitsT>
void GlobalStoichiometricInteraction<traitsT>::calc_force(system_type& sys) const noexcept
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();

    // initialization of each buffering container.
    std::fill(pots_sum_a_.begin(), pots_sum_a_.end(), 0.0);
    std::fill(pots_sum_b_.begin(), pots_sum_b_.end(), 0.0);
    std::fill(derivs_sum_a_.begin(), derivs_sum_a_.end(),
              math::make_coordinate<coordinate_type>(0.0, 0.0, 0.0));
    std::fill(derivs_sum_b_.begin(), derivs_sum_b_.end(),
              math::make_coordinate<coordinate_type>(0.0, 0.0, 0.0));

    // make index pair list to specify the range of buffer for specific a and
    // pre calculation for pots_sum and derivs_sum for each particle.
    const auto        participants_a       = potential_.participants_a();
    const std::size_t participants_a_num   = potential_.participants_a_num();

    for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
    {
        const index_type i = participants_a[idx_a];

        partner_buffer_tuple_type& partner_buffer_range = partner_buffer_ranges_[idx_a];
        const auto& partner     = std::get<0>(partner_buffer_range);
        const auto& pot_range   = std::get<1>(partner_buffer_range);
        const auto& deriv_range = std::get<2>(partner_buffer_range);

        real_type&       pots_sum_a   = pots_sum_a_  [idx_a];
        coordinate_type& derivs_sum_a = derivs_sum_a_[idx_a];
        for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
        {
            const index_type      j     = partner[ptnr_idx].index;
            const coordinate_type rij   = sys.adjust_direction(sys.position(i), sys.position(j));
            const real_type       l2    = math::length_sq(rij); // |rij|^2
            const real_type       rl    = math::rsqrt(l2);      // 1 / |rij|
            const real_type       l     = l2 * rl;
            const coordinate_type deriv = potential_.derivative(l) * rl * rij;
            const real_type       pot   = potential_.potential(l);

            pot_range  [ptnr_idx] =  pot;
            deriv_range[ptnr_idx] =  deriv;
            pots_sum_a            += pot;
            derivs_sum_a          += deriv;

            const index_type idx_b  =  idx_buffer_map_[j];
            pots_sum_b_      [idx_b] += pot;
            derivs_sum_b_    [idx_b] -= deriv;
        }
    }

    // force calculation
    for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
    {
        const index_type i = participants_a[idx_a];

        const auto&     partner_buffer_tuple  = partner_buffer_ranges_[idx_a];
        const auto&     partner               = std::get<0>(partner_buffer_tuple);
        const auto&     potential_buff_range  = std::get<1>(partner_buffer_tuple);
        const auto&     derivative_buff_range = std::get<2>(partner_buffer_tuple);
        const real_type pots_sum_a             = pots_sum_a_[idx_a];

        if(pots_sum_a <= 1.0)
        {
            for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
            {
                const index_type       j             = partner[ptnr_idx].index;
                const index_type       idx_b         = idx_buffer_map_[j];
                const real_type        pots_sum_b    = pots_sum_b_[idx_b];
                const coordinate_type& deriv_buff_ab = derivative_buff_range[ptnr_idx];

                if(pots_sum_b <= 1.0)
                {
                    const coordinate_type force = epsilon_ * deriv_buff_ab;
                    sys.force(i) -= force;
                    sys.force(j) += force;
                }
                else // pots_sum_a < pots_sum_b
                {
                    const real_type        pot_buff_ab    = potential_buff_range[ptnr_idx];
                    const real_type        pots_sum_b     = pots_sum_b_         [idx_b];
                    const coordinate_type& derivs_sum_b   = derivs_sum_b_       [idx_b];
                    const real_type        inv_pots_sum_b = 1.0 / pots_sum_b;
                    const real_type        ep_pots_sum_b  = epsilon_ * inv_pots_sum_b;
                    const real_type        ep_pots_sum_b2 = ep_pots_sum_b * inv_pots_sum_b;

                    sys.force(i) -= ep_pots_sum_b * deriv_buff_ab;
                    sys.force(j) +=
                        ep_pots_sum_b2 * (pot_buff_ab * derivs_sum_b + pots_sum_b * deriv_buff_ab);

                    // the force to not a or b particle from a-b pair
                    const auto& j_partner = partition_.partners(j);
                    for(const auto& j_ptnr : j_partner)
                    {
                        const index_type k     = j_ptnr.index;
                        const index_type idx_c = idx_buffer_map_[k];
                        const auto&      partner_buffer_tuple_c = partner_buffer_ranges_[idx_c];
                        const auto&      deriv_buff_range_c     = std::get<2>(partner_buffer_tuple_c);

                        sys.force(k) +=
                            ep_pots_sum_b2 * pot_buff_ab * deriv_buff_range_c[idx_b];
                    }
                }
            }
        }
        else // 1.0 < pots_sum_a case
        {
            const real_type inv_pots_sum_a      = 1.0 / pots_sum_a;
            const real_type ep_pots_sum_a       = epsilon_ * inv_pots_sum_a;
            const real_type ep_pots_sum_a2      = ep_pots_sum_a * inv_pots_sum_a;
            const coordinate_type& derivs_sum_a = derivs_sum_a_[idx_a];
            const coordinate_type  derivs_sum_a_pots_sum_a = inv_pots_sum_a * derivs_sum_a;

            for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
            {
                const index_type j     = partner[ptnr_idx].index;
                const index_type idx_b = idx_buffer_map_[j];
                const real_type        pot_buff_ab   = potential_buff_range [ptnr_idx];
                const coordinate_type& deriv_buff_ab = derivative_buff_range[ptnr_idx];
                const real_type        pots_sum_b    = pots_sum_b_          [idx_b];
                const coordinate_type& derivs_sum_b  = derivs_sum_b_        [idx_b];

                if(pots_sum_b < pots_sum_a)
                {
                    sys.force(i) +=
                        ep_pots_sum_a2 * (derivs_sum_a * pot_buff_ab - pots_sum_a * deriv_buff_ab);
                    sys.force(j) += ep_pots_sum_a  * deriv_buff_ab;

                    // the force to not a or b particle from a-b pair
                    for(const auto& i_ptnr : partner)
                    {
                        const index_type k     = i_ptnr.index;
                        const index_type idx_d = idx_buffer_map_[k];

                        sys.force(k) -=
                            ep_pots_sum_a2 * pot_buff_ab * derivative_buff_range[idx_d];
                    }
                }
                else // pots_sum_a <= pots_sum_b
                {
                    const real_type inv_pots_sum_b  = 1.0 / pots_sum_b;
                    const real_type ep_pots_sum_b   = epsilon_ * inv_pots_sum_b;
                    const real_type ep_pots_sum_b2  = ep_pots_sum_b * inv_pots_sum_b;

                    sys.force(i) -= ep_pots_sum_b * deriv_buff_ab;
                    sys.force(j) +=
                        ep_pots_sum_b2 * (pot_buff_ab * derivs_sum_b + pots_sum_b * deriv_buff_ab);

                    // the force to not a or b particle from a-b pair
                    const auto& j_partner = partition_.partners(j);
                    for(const auto& j_ptnr : j_partner)
                    {
                        const index_type k     = j_ptnr.index;
                        const index_type idx_c = idx_buffer_map_[k];
                        const auto&      partner_buffer_tuple_c = partner_buffer_ranges_[idx_c];
                        const auto&      deriv_buff_range_c     = std::get<2>(partner_buffer_tuple_c);

                        sys.force(k) +=
                            ep_pots_sum_b2 * pot_buff_ab * deriv_buff_range_c[idx_b];
                    }
                }
            }
        }
    }
    return;
}

template<typename traitsT>
typename GlobalStoichiometricInteraction<traitsT>::real_type
GlobalStoichiometricInteraction<traitsT>::calc_energy(const system_type& sys) const noexcept
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();

    // initialization of each buffering container.
    std::fill(pots_sum_a_.begin(), pots_sum_a_.end(), 0.0);
    std::fill(pots_sum_b_.begin(), pots_sum_b_.end(), 0.0);

    // make index pair list to specify the range of buffer for specific a and
    // pre calculation for pots_sum and derivs_sum for each particle.
    const auto        participants_a     = potential_.participants_a();
    const std::size_t participants_a_num = potential_.participants_a_num();

    for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
    {
        const index_type   i = participants_a[idx_a];

        partner_buffer_tuple_type& partner_buffer_range = partner_buffer_ranges_[idx_a];
        const auto& partner   = std::get<0>(partner_buffer_range);
        const auto& pot_range = std::get<1>(partner_buffer_range);

        real_type& pots_sum_a = pots_sum_a_[idx_a];

        for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
        {
            const index_type j     = partner[ptnr_idx].index;
            const index_type idx_b = idx_buffer_map_[j];
            const coordinate_type rij =
                sys.adjust_direction(sys.position(i), sys.position(j));
            const real_type l2  = math::length_sq(rij); // |rij|^2
            const real_type rl  = math::rsqrt(l2);      // 1 / |rij|
            const real_type l   = l2 * rl;              // |rij|
            const real_type pot = potential_.potential(l);

            pot_range[ptnr_idx] =  pot;
            pots_sum_a           += pot;
            pots_sum_b_  [idx_b] += pot;
        }
    }

    // potential calculation
    real_type retval(0.0);
    for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
    {
        const partner_buffer_tuple_type&   partner_buffer_range = partner_buffer_ranges_[idx_a];
        const partner_range_type&          partner              = std::get<0>(partner_buffer_range);
        const potential_buffer_range_type& pots_buff_a          = std::get<1>(partner_buffer_range);
        for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
        {
            const index_type j            = partner[ptnr_idx].index;
            const index_type idx_b        = idx_buffer_map_[j];
            const real_type  pot_buff_ab  = pots_buff_a[ptnr_idx];
            const real_type  max_pots_sum  = std::max(pots_sum_a_[idx_a], pots_sum_b_[idx_b]);
            if(max_pots_sum <= 1.0)
            {
                retval += pot_buff_ab;
            }
            else
            {
                retval += pot_buff_ab / max_pots_sum;
            }
        }
    }
    retval *= -epsilon_;

    return retval;
}

template<typename traitsT>
typename GlobalStoichiometricInteraction<traitsT>::real_type
GlobalStoichiometricInteraction<traitsT>::calc_force_and_energy(system_type& sys) const noexcept
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();

    // initialization of each buffering container.
    std::fill(pots_sum_a_.begin(), pots_sum_a_.end(), 0.0);
    std::fill(pots_sum_b_.begin(), pots_sum_b_.end(), 0.0);
    std::fill(derivs_sum_a_.begin(), derivs_sum_a_.end(),
              math::make_coordinate<coordinate_type>(0.0, 0.0, 0.0));
    std::fill(derivs_sum_b_.begin(), derivs_sum_b_.end(),
              math::make_coordinate<coordinate_type>(0.0, 0.0, 0.0));

    // make index pair list to specify the range of buffer for specific a and
    // pre calculation for pots_sum and derivs_sumb for each particle.
    const auto        participants_a     = potential_.participants_a();
    const std::size_t participants_a_num = potential_.participants_a_num();

    for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
    {
        const index_type i = participants_a[idx_a];

        partner_buffer_tuple_type& partner_buffer_range = partner_buffer_ranges_[idx_a];
        const auto& partner     = std::get<0>(partner_buffer_range);
        const auto& pot_range   = std::get<1>(partner_buffer_range);
        const auto& deriv_range = std::get<2>(partner_buffer_range);

        real_type&       pots_sum_a   = pots_sum_a_  [idx_a];
        coordinate_type& derivs_sum_a = derivs_sum_a_[idx_a];
        for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
        {
            const index_type      j     = partner[ptnr_idx].index;
            const coordinate_type rij   = sys.adjust_direction(sys.position(i), sys.position(j));
            const real_type       l2    = math::length_sq(rij); // |rij|^2
            const real_type       rl    = math::rsqrt(l2);      // 1 / |rij|
            const real_type       l     = l2 * rl;
            const coordinate_type deriv = potential_.derivative(l) * rl * rij;
            const real_type       pot   = potential_.potential(l);

            pot_range    [ptnr_idx] =  pot;
            deriv_range  [ptnr_idx] =  deriv;
            pots_sum_a              += pot;
            derivs_sum_a            += deriv;

            const index_type idx_b =  idx_buffer_map_[j];
            pots_sum_b_  [idx_b]   += pot;
            derivs_sum_b_[idx_b]   -= deriv;
        }
    }

    // force calculation
    for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
    {
        const index_type i = participants_a[idx_a];

        const auto&     partner_buffer_tuple  = partner_buffer_ranges_[idx_a];
        const auto&     partner               = std::get<0>(partner_buffer_tuple);
        const auto&     potential_buff_range  = std::get<1>(partner_buffer_tuple);
        const auto&     derivative_buff_range = std::get<2>(partner_buffer_tuple);
        const real_type pots_sum_a            = pots_sum_a_[idx_a];

        if(pots_sum_a <= 1.0)
        {
            for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
            {
                const index_type       j             = partner[ptnr_idx].index;
                const index_type       idx_b         = idx_buffer_map_[j];
                const real_type        pots_sum_b    = pots_sum_b_  [idx_b];
                const coordinate_type& derivs_sum_b  = derivs_sum_b_[idx_b];
                const coordinate_type& deriv_buff_ab = derivative_buff_range[ptnr_idx];

                if(pots_sum_b <= 1.0)
                {
                    const coordinate_type force = epsilon_ * deriv_buff_ab;
                    sys.force(i) -= force;
                    sys.force(j) += force;
                }
                else // pots_sum_a < pots_sum_b case
                {
                    const real_type        pot_buff_ab    = potential_buff_range [ptnr_idx];
                    const coordinate_type& deriv_buff_ab  = derivative_buff_range[ptnr_idx];
                    const real_type        inv_pots_sum_b = 1.0 / pots_sum_b;
                    const real_type        ep_pots_sum_b  = epsilon_ * inv_pots_sum_b;
                    const real_type        ep_pots_sum_b2 = ep_pots_sum_b * inv_pots_sum_b;

                    sys.force(i) -= ep_pots_sum_b * deriv_buff_ab;
                    sys.force(j) +=
                        ep_pots_sum_b2 * (pot_buff_ab * derivs_sum_b + pots_sum_b * deriv_buff_ab);

                    // the force to not a or b particle from a-b pair
                    const auto& j_partner = partition_.partners(j);
                    for(const auto& j_ptnr : j_partner)
                    {
                        const index_type k     = j_ptnr.index;
                        const index_type idx_c = idx_buffer_map_[k];
                        const auto&      partner_buffer_tuple_c = partner_buffer_ranges_[idx_c];
                        const auto&      deriv_buff_range_c     = std::get<2>(partner_buffer_tuple_c);

                        sys.force(k) +=
                            ep_pots_sum_b2 * pot_buff_ab * deriv_buff_range_c[idx_b];
                    }
                }
            }
        }
        else // 1.0 < pots_sum_a case
        {
            const auto&            potential_buff_range = std::get<1>(partner_buffer_tuple);
            const real_type        inv_pots_sum_a       = 1.0 / pots_sum_a;
            const real_type        ep_pots_sum_a        = epsilon_ * inv_pots_sum_a;
            const real_type        ep_pots_sum_a2       = ep_pots_sum_a * inv_pots_sum_a;
            const coordinate_type& derivs_sum_a         = derivs_sum_a_[idx_a];

            for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
            {
                const index_type j     = partner[ptnr_idx].index;
                const index_type idx_b = idx_buffer_map_[j];
                const real_type        pot_buff_ab   = potential_buff_range [ptnr_idx];
                const coordinate_type& deriv_buff_ab = derivative_buff_range[ptnr_idx];
                const real_type        pots_sum_b    = pots_sum_b_          [idx_b];
                const coordinate_type& derivs_sum_b  = derivs_sum_b_        [idx_b];

                if(pots_sum_b < pots_sum_a)
                {
                    sys.force(i) += 
                        ep_pots_sum_a2 * (derivs_sum_a * pot_buff_ab - pots_sum_a * deriv_buff_ab);
                    sys.force(j) += ep_pots_sum_a * deriv_buff_ab;

                    // the force to not a or b particle from a-b pair
                    for(const auto& i_ptnr : partner)
                    {
                        const index_type k     = i_ptnr.index;
                        const index_type idx_d = idx_buffer_map_[k];

                        sys.force(k) -=
                            ep_pots_sum_a2 * pot_buff_ab * derivative_buff_range[idx_d];
                    }
                }
                else // pots_sum_a_ <= pots_sum_b_
                {
                    const real_type inv_pots_sum_b = 1.0 / pots_sum_b;
                    const real_type ep_pots_sum_b  = epsilon_ * inv_pots_sum_b;
                    const real_type ep_pots_sum_b2 = ep_pots_sum_b * inv_pots_sum_b;

                    sys.force(i) += ep_pots_sum_b * deriv_buff_ab;
                    sys.force(j) -=
                        ep_pots_sum_b * (pot_buff_ab * derivs_sum_b + pots_sum_b * deriv_buff_ab);

                    // the force to not a or b particle from a-b pair
                    const auto& j_partner = partition_.partners(j);
                    for(const auto& j_ptnr : j_partner)
                    {
                        const index_type k     = j_ptnr.index;
                        const index_type idx_c = idx_buffer_map_[k];
                        const auto&      partner_buffer_tuple_c = partner_buffer_ranges_[idx_c];
                        const auto&      deriv_buff_range_c     = std::get<2>(partner_buffer_tuple_c);

                        sys.force(k) +=
                            ep_pots_sum_b2 * pot_buff_ab * deriv_buff_range_c[idx_b];
                    }
                }
            }
        }
    }

    // potential calculation
    real_type retval(0.0);
    for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a) 
    {
        partner_buffer_tuple_type&   partner_buffer_range = partner_buffer_ranges_[idx_a];
        partner_range_type&          partner              = std::get<0>(partner_buffer_range);
        potential_buffer_range_type& pots_buff_a          = std::get<1>(partner_buffer_range);
        for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
        {
            const index_type j            = partner[ptnr_idx].index;
            const index_type idx_b        = idx_buffer_map_[j];
            const real_type  pots_buff_ab = pots_buff_a[ptnr_idx];
            const real_type  max_pots_sum = std::max(pots_sum_a_[idx_a], pots_sum_b_[idx_b]);
            if(max_pots_sum <= 1.0)
            {
                retval += pots_buff_ab;
            }
            else
            {
                retval += pots_buff_ab / max_pots_sum;
            }
        }
    }
    retval *= -epsilon_;

    return retval;
}
} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize major use-cases
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

namespace mjolnir
{
extern template class GlobalStoichiometricInteraction<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class GlobalStoichiometricInteraction<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class GlobalStoichiometricInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class GlobalStoichiometricInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_FORCEFIELD_STOICHIOMETRIC_GLOBAL_STOICHIOMETRIC_INTERACTION_HPP */
