#ifndef MJOLNIR_FORCEFIELD_STOICHIOMETRIC_GLOBAL_STOICHIOMETRIC_INTERACTION_HPP
#define MJOLNIR_FORCEFIELD_STOICHIOMETRIC_GLOBAL_STOICHIOMETRIC_INTERACTION_HPP
#include <mjolnir/forcefield/global/ParameterList.hpp>
#include <mjolnir/forcefield/stoichiometric/GlobalStoichiometricInteractionPotential.hpp>
#include <mjolnir/forcefield/stoichiometric/StoichiometricInteractionRule.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/GlobalInteractionBase.hpp>
#include <mjolnir/core/SpatialPartitionBase.hpp>
#include <mjolnir/util/logger.hpp>
#include <memory>

namespace mjolnir
{

// The interaction which conserve stoichiometry.

template<typename traitsT, typename potentialT>
class GlobalStoichiometricInteraction final : public GlobalInteractionBase<traitsT>
{
  public:

    using traits_type            = traitsT;
    using potential_type         = potentialT;
    using base_type              = GlobalInteractionBase<traitsT>;
    using index_type             = std::size_t;
    using real_type              = typename base_type::real_type;
    using coordinate_type        = typename base_type::coordinate_type;
    using system_type            = typename base_type::system_type;
    using topology_type          = typename base_type::topology_type;
    using partition_type         = SpatialPartition<traitsT, potential_type>;
    using parameter_list_type    = ParameterList<traits_type, potential_type>;
    using potential_buffer_type  = std::vector<real_type>;
    using derivative_buffer_type = std::vector<coordinate_type>;

    using potential_buffer_iterator    = typename potential_buffer_type::iterator;
    using derivative_buffer_iterator   = typename derivative_buffer_type::iterator;
    using partner_range_type           = typename partition_type::range_type;
    using partner_buffer_type          = std::vector<partner_range_type>;
    using potential_range_type         = range<potential_buffer_iterator>;
    using derivative_range_type        = range<derivative_buffer_iterator>;
    using potential_range_buffer_type  = std::vector<potential_range_type>;
    using derivative_range_buffer_type = std::vector<derivative_range_type>;

  public:
    GlobalStoichiometricInteraction(potential_type&& pot, parameter_list_type&& para,
        partition_type&& part,
        real_type epsilon, std::size_t coefa, std::size_t coefb)
        : potential_(std::move(pot)), parameters_(std::move(para)),
          partition_(std::move(part)),
          epsilon_(epsilon), coefa_(coefa), coefb_(coefb), inv_coefa_(1./coefa), inv_coefb_(1./coefb)
    {
        const std::size_t participants_a_num = parameters_.participants_a_num();
        const std::size_t participants_b_num = parameters_.participants_b_num();
        const std::size_t participants_ab    = participants_a_num * participants_b_num;
        derivs_buff_    .resize(participants_ab);
        potentials_buff_.resize(participants_ab);
        partner_buffer_a_    .resize(participants_a_num);
        pot_buffer_range_a_  .resize(participants_a_num);
        deriv_buffer_range_a_.resize(participants_a_num);

        pots_sum_a_.resize(participants_a_num);
        pots_sum_b_.resize(participants_b_num);
        pots_partial_sum_a_.resize(participants_a_num);
        pots_partial_sum_b_.resize(participants_b_num);
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

        this->potential_ .initialize(sys);
        this->parameters_.initialize(sys, topol, potential_);
        this->partition_ .initialize(sys, parameters_.cref());

        this->update_buffer_range();

        // When we calculate potential or derivative for specific idx, we buffering the
        // value in intermediate value container. So, we have to know which index in
        // the system correspond to which index of buffering container. These for loop
        // is to make that container.
        idx_buffer_map_.resize(sys.size());
        const auto participants_a = this->parameters_.participants_a();
        const auto participants_b = this->parameters_.participants_b();
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
        this->potential_ .update(sys);
        this->parameters_.update(sys, topol, potential_);
        this->partition_ .initialize(sys, this->parameters_cref());

        this->update_buffer_range();
    }

    void reduce_margin(const real_type dmargin, const system_type& sys) override
    {
        if(partition_.reduce_margin(dmargin, sys, this->parametes_.cref()))
        {
            this->update_buffer_range();
        }
        return;
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        partition_.scale_margin(scale, sys, this->parameters_.cref());
        return;
    }

    void      calc_force (system_type&)           const noexcept override;
    void      calc_force_and_virial(system_type&) const noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FATAL("GlobalStoichiometricInteraction:"
                " virial calculation is not supported");
    }
    real_type calc_force_and_energy(system_type&) const noexcept override;
    real_type calc_force_virial_energy(system_type&) const noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FATAL("GlobalStoichiometricInteraction:"
                " virial calculation is not supported");
    }
    real_type calc_energy(const system_type&)     const noexcept override;

    std::string name() const override {return "Stoichiometric";}

    potential_type const& potential() const noexcept {return potential_;}
    potential_type &      potential()       noexcept {return potential_;}

    base_type* clone() const override
    {
        return new GlobalStoichiometricInteraction(
                potential_type(potential_), parameter_list_type(parameters_),
                epsilon_, coefa_, coefb_);
    }

  private:
    void update_buffer_range()
    {
        MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
        const std::vector<index_type> participants_a
            = parameters_.participants_a();
        const std::size_t participants_a_num = parameters_.participants_a_num();
        auto pot_buff_iter   = potentials_buff_.begin();
        auto deriv_buff_iter = derivs_buff_    .begin();

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

            partner_buffer_a_    [idx_a] = partner;
            pot_buffer_range_a_  [idx_a] = pot_range;
            deriv_buffer_range_a_[idx_a] = deriv_range;
        }
    }

  private:

    potential_type      potential_;
    parameter_list_type parameters_;
    partition_type      partition_;

    real_type      epsilon_;
    std::size_t    coefa_;
    std::size_t    coefb_;
    real_type      inv_coefa_;
    real_type      inv_coefb_;

    // -----------------------------------------------------------------------
    // Variable for mapping system index to buffering index
    std::vector<std::size_t> idx_buffer_map_;

    // Variables for buffering intermediate value.
    // These value can change in calc_force and calc_energy function.
    mutable potential_buffer_type  potentials_buff_;
    mutable derivative_buffer_type derivs_buff_;
    mutable potential_buffer_type  pots_sum_a_;
    mutable potential_buffer_type  pots_sum_b_;
    mutable potential_buffer_type  pots_partial_sum_a_;
    mutable potential_buffer_type  pots_partial_sum_b_;

    // Variables to specify the range of specific a partners in buffer.
    mutable partner_buffer_type          partner_buffer_a_;
    mutable potential_range_buffer_type  pot_buffer_range_a_;
    mutable derivative_range_buffer_type deriv_buffer_range_a_;

    // sum of derivation of potential function for specific first particle.
    mutable derivative_buffer_type     derivs_sum_a_;

    // sum of derivation of potential function for specific second particle.
    mutable derivative_buffer_type     derivs_sum_b_;
};

template<typename traitsT, typename potT>
void GlobalStoichiometricInteraction<traitsT, potT>::calc_force(system_type& sys) const noexcept
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();

    // initialization of each buffering container.
    std::fill(pots_sum_a_.begin(),         pots_sum_a_.end(),         0.0);
    std::fill(pots_partial_sum_a_.begin(), pots_partial_sum_a_.end(), 0.0);
    std::fill(pots_sum_b_.begin(),         pots_sum_b_.end(),         0.0);
    std::fill(pots_partial_sum_b_.begin(), pots_partial_sum_b_.end(), 0.0);
    std::fill(derivs_sum_a_.begin(), derivs_sum_a_.end(),
              math::make_coordinate<coordinate_type>(0.0, 0.0, 0.0));
    std::fill(derivs_sum_b_.begin(), derivs_sum_b_.end(),
              math::make_coordinate<coordinate_type>(0.0, 0.0, 0.0));

    // make index pair list to specify the range of buffer for specific a and
    // pre calculation for pots_sum and derivs_sum for each particle.
    const auto&       participants_a     = parameters_.participants_a();
    const std::size_t participants_a_num = parameters_.participants_a_num();

    for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
    {
        const index_type i = participants_a[idx_a];

        const auto& partner     = partner_buffer_a_    [idx_a];
        const auto& pot_range   = pot_buffer_range_a_  [idx_a];
        const auto& deriv_range = deriv_buffer_range_a_[idx_a];

        real_type&       pots_sum_a   = pots_sum_a_  [idx_a];
        coordinate_type& derivs_sum_a = derivs_sum_a_[idx_a];
        for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
        {
            const auto&           ptnr  = partner[ptnr_idx];
            const index_type      j     = ptnr.index;
            const auto&           para  = ptnr.parameter();

            const coordinate_type rij   = sys.adjust_direction(sys.position(i), sys.position(j));
            const real_type       l2    = math::length_sq(rij); // |rij|^2
            const real_type       rl    = math::rsqrt(l2);      // 1 / |rij|
            const real_type       l     = l2 * rl;
            const coordinate_type deriv = potential_.derivative(l, para) * rl * rij;
            const real_type       pot   = potential_.potential(l, para);

            pot_range  [ptnr_idx] =  pot;
            deriv_range[ptnr_idx] =  deriv;
            pots_sum_a            += pot;
            derivs_sum_a          += deriv;

            const index_type  idx_b  =  idx_buffer_map_[j];
            pots_sum_b_      [idx_b] += pot;
            derivs_sum_b_    [idx_b] -= deriv;
        }
    }

    for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
    {
        const auto&     partner            = partner_buffer_a_  [idx_a];
        const auto&     pot_buff_range     = pot_buffer_range_a_[idx_a];
        const real_type pots_sum_a         = pots_sum_a_        [idx_a];
        const real_type pots_sum_a_coefb   = pots_sum_a*inv_coefb_;
        real_type&      pots_partial_sum_a = pots_partial_sum_a_[idx_a];

        for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
        {
            const index_type j                = partner[ptnr_idx].index;
            const index_type idx_b            = idx_buffer_map_[j];
            const real_type  pots_sum_b       = pots_sum_b_[idx_b];
            const real_type  pots_sum_b_coefa = pots_sum_b*inv_coefa_;

            if(std::max(pots_sum_a_coefb, pots_sum_b_coefa) <= 1.0){ continue; }

            const real_type  pot_buff_ab = pot_buff_range[ptnr_idx];
            if(pots_sum_a_coefb < pots_sum_b_coefa)
            {
                pots_partial_sum_b_[idx_b] += pot_buff_ab;
            }
            else // pots_sum_a / coefb >= pots_sum_b / coefa
            {
                pots_partial_sum_a += pot_buff_ab;
            }
        }
    }

    // force calculation
    for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
    {
        const index_type i = participants_a[idx_a];

        const auto& partner               = partner_buffer_a_    [idx_a];
        const auto& potential_buff_range  = pot_buffer_range_a_  [idx_a];
        const auto& derivative_buff_range = deriv_buffer_range_a_[idx_a];

        const real_type        pots_sum_a         = pots_sum_a_        [idx_a];
        const coordinate_type& derivs_sum_a       = derivs_sum_a_      [idx_a];
        const real_type        pots_partial_sum_a = pots_partial_sum_a_[idx_a];
        const real_type        inv_pots_sum_a     = 1.0 / pots_sum_a;
        const real_type        coefb_pots_sum_a   = coefb_ * inv_pots_sum_a;
        const real_type        ep_coefb_pots_sum_a  = epsilon_ * coefb_pots_sum_a;
        const real_type        ep_coefb_pots_sum_a2 = ep_coefb_pots_sum_a * inv_pots_sum_a;

        for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
        {
            const index_type       j             = partner[ptnr_idx].index;
            const index_type       idx_b         = idx_buffer_map_[j];
            const real_type        pot_buff_ab   = potential_buff_range [ptnr_idx];
            const coordinate_type& deriv_buff_ab = derivative_buff_range[ptnr_idx];
            const real_type        pots_sum_b    = pots_sum_b_[idx_b];
            const real_type        inv_pots_sum_b   = 1.0 / pots_sum_b;
            const real_type        coefa_pots_sum_b = coefa_ * inv_pots_sum_b;
            const coordinate_type& derivs_sum_b     = derivs_sum_b_[idx_b];

            if(1.0 <= coefb_pots_sum_a) // pots_sum_a / coefb <= 1.0 case
            {
                if(1.0 <= coefa_pots_sum_b) // pots_sum_b / coefa <= 1.0 case
                {
                    // force from distance part derivation
                    const coordinate_type force = epsilon_ * deriv_buff_ab;
                    sys.force(i) -= force;
                    sys.force(j) += force;
                }
                else // pots_sum_a / coefb < pots_sum_b / coefa
                {
                    const real_type ep_coefa_pots_sum_b  = epsilon_ * coefa_pots_sum_b;
                    const real_type ep_coefa_pots_sum_b2 = ep_coefa_pots_sum_b * inv_pots_sum_b;

                    // force from distance part derivation
                    const coordinate_type dforce =
                        ep_coefa_pots_sum_b * deriv_buff_ab;
                    sys.force(i) -= dforce;
                    sys.force(j) += dforce;

                    // force from valence part derivation
                    const coordinate_type vforce =
                        ep_coefa_pots_sum_b2 * derivs_sum_b * pot_buff_ab;
                    sys.force(j) += vforce;
                }
            }
            else // 1.0 < pots_sum_a / coefb case
            {
                // force from depends on valece part
                if(coefb_pots_sum_a <= coefa_pots_sum_b) // pots_sum_b / coefa <= pots_sum_a / coefb
                {
                    // force from distance part derivation
                    sys.force(i) -= ep_coefb_pots_sum_a * deriv_buff_ab;
                    sys.force(j) += ep_coefb_pots_sum_a * deriv_buff_ab;

                    // force from valence part derivation
                    sys.force(i) += ep_coefb_pots_sum_a2 * derivs_sum_a * pot_buff_ab;
                }
                else // pots_sum_a / coefb < pots_sum_b / coefa
                {
                    const real_type ep_coefa_pots_sum_b  = epsilon_* coefa_pots_sum_b;
                    const real_type ep_coefa_pots_sum_b2 = ep_coefa_pots_sum_b * inv_pots_sum_b;

                    // force from distance part derivation
                    sys.force(i) -= ep_coefa_pots_sum_b * deriv_buff_ab;
                    sys.force(j) += ep_coefa_pots_sum_b * deriv_buff_ab;
                    // force from valence part derivation
                    sys.force(j) += ep_coefa_pots_sum_b2 * derivs_sum_b * pot_buff_ab;
                }
                // force from depends on valence part
                sys.force(j) -= ep_coefb_pots_sum_a2 * deriv_buff_ab * pots_partial_sum_a;
            }
        }

        // force from depends on valence part
        for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
        {
            const index_type k          = partner[ptnr_idx].index;
            const index_type idx_d      = idx_buffer_map_[k];
            const real_type  pots_sum_d = pots_sum_b_[idx_d];
            const real_type  inv_pots_sum_d   = 1.0 / pots_sum_d;
            const real_type  coefa_pots_sum_d = coefa_ * inv_pots_sum_d;

            if(coefa_pots_sum_d < 1.0)
            {
                const real_type        ep_coefa_pots_sum_d  = epsilon_ * coefa_pots_sum_d;
                const real_type        ep_coefa_pots_sum_d2 = ep_coefa_pots_sum_d * inv_pots_sum_d;
                const coordinate_type& deriv_buff_ad  = derivative_buff_range[ptnr_idx];

                sys.force(i) += ep_coefa_pots_sum_d2 * deriv_buff_ad * pots_partial_sum_b_[idx_d];
            }
        }
    }
    return;
}

template<typename traitsT, typename potT>
typename GlobalStoichiometricInteraction<traitsT, potT>::real_type
GlobalStoichiometricInteraction<traitsT, potT>::calc_energy(const system_type& sys) const noexcept
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();

    // initialization of each buffering container.
    std::fill(pots_sum_a_.begin(), pots_sum_a_.end(), 0.0);
    std::fill(pots_sum_b_.begin(), pots_sum_b_.end(), 0.0);

    // make index pair list to specify the range of buffer for specific a and
    // pre calculation for pots_sum and derivs_sum for each particle.
    const auto        participants_a     = parameters_.participants_a();
    const std::size_t participants_a_num = parameters_.participants_a_num();

    for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
    {
        const index_type i  = participants_a[idx_a];

        const auto& partner   = partner_buffer_a_  [idx_a];
        const auto& pot_range = pot_buffer_range_a_[idx_a];

        real_type& pots_sum_a = pots_sum_a_[idx_a];

        for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
        {
            const auto&      ptnr  = partner[ptnr_idx].index;
            const index_type j     = ptnr.index;
            const auto&      para  = ptnr.parameter();

            const index_type idx_b = idx_buffer_map_[j];
            const coordinate_type rij =
                sys.adjust_direction(sys.position(i), sys.position(j));
            const real_type l2  = math::length_sq(rij); // |rij|^2
            const real_type rl  = math::rsqrt(l2);      // 1 / |rij|
            const real_type l   = l2 * rl;              // |rij|
            const real_type pot = potential_.potential(l, para);

            pot_range [ptnr_idx] =  pot;
            pots_sum_a           += pot;
            pots_sum_b_  [idx_b] += pot;
        }
    }

    // potential calculation
    real_type retval(0.0);
    for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
    {
        const partner_range_type&   partner     = partner_buffer_a_  [idx_a];
        const potential_range_type& pots_buff_a = pot_buffer_range_a_[idx_a];
        for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
        {
            const index_type j            = partner[ptnr_idx].index;
            const index_type idx_b        = idx_buffer_map_[j];
            const real_type  pot_buff_ab  = pots_buff_a[ptnr_idx];
            const real_type  max_normalized_contact_order =
                std::max(pots_sum_a_[idx_a]*inv_coefb_, pots_sum_b_[idx_b]*inv_coefa_);
            if(max_normalized_contact_order <= 1.0)
            {
                retval += pot_buff_ab;
            }
            else
            {
                retval += pot_buff_ab / max_normalized_contact_order;
            }
        }
    }
    retval *= -epsilon_;

    return retval;
}

template<typename traitsT, typename potT>
typename GlobalStoichiometricInteraction<traitsT, potT>::real_type
GlobalStoichiometricInteraction<traitsT, potT>::calc_force_and_energy(system_type& sys) const noexcept
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();

    // initialization of each buffering container.
    std::fill(pots_sum_a_.begin(),         pots_sum_a_.end(),         0.0);
    std::fill(pots_partial_sum_a_.begin(), pots_partial_sum_a_.end(), 0.0);
    std::fill(pots_sum_b_.begin(),         pots_sum_b_.end(),         0.0);
    std::fill(pots_partial_sum_b_.begin(), pots_partial_sum_b_.end(), 0.0);
    std::fill(derivs_sum_a_.begin(), derivs_sum_a_.end(),
              math::make_coordinate<coordinate_type>(0.0, 0.0, 0.0));
    std::fill(derivs_sum_b_.begin(), derivs_sum_b_.end(),
              math::make_coordinate<coordinate_type>(0.0, 0.0, 0.0));

    // make index pair list to specify the range of buffer for specific a and
    // pre calculation for pots_sum and derivs_sumb for each particle.
    const auto        participants_a     = parameters_.participants_a();
    const std::size_t participants_a_num = parameters_.participants_a_num();

    for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
    {
        const index_type i = participants_a[idx_a];

        const auto& partner     = partner_buffer_a_    [idx_a];
        const auto& pot_range   = pot_buffer_range_a_  [idx_a];
        const auto& deriv_range = deriv_buffer_range_a_[idx_a];

        real_type&       pots_sum_a   = pots_sum_a_  [idx_a];
        coordinate_type& derivs_sum_a = derivs_sum_a_[idx_a];
        for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
        {
            const auto&           ptnr  = partner[ptnr_idx];
            const index_type      j     = ptnr.index;
            const auto&           para  = ptnr.parameter();;;;

            const coordinate_type rij   = sys.adjust_direction(sys.position(i), sys.position(j));
            const real_type       l2    = math::length_sq(rij); // |rij|^2
            const real_type       rl    = math::rsqrt(l2);      // 1 / |rij|
            const real_type       l     = l2 * rl;
            const coordinate_type deriv = potential_.derivative(l, para) * rl * rij;
            const real_type       pot   = potential_.potential(l, para);

            pot_range    [ptnr_idx] =  pot;
            deriv_range  [ptnr_idx] =  deriv;
            pots_sum_a              += pot;
            derivs_sum_a            += deriv;

            const index_type idx_b =  idx_buffer_map_[j];
            pots_sum_b_  [idx_b]   += pot;
            derivs_sum_b_[idx_b]   -= deriv;
        }
    }

    for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
    {
        const auto&      partner            = partner_buffer_a_  [idx_a];
        const auto&      pot_buff_range     = pot_buffer_range_a_[idx_a];
        const real_type  pots_sum_a         = pots_sum_a_[idx_a];
        const real_type  pots_sum_a_coefb   = pots_sum_a*inv_coefb_;
        real_type&       pots_partial_sum_a = pots_partial_sum_a_[idx_a];

        for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
        {
            const index_type j                = partner[ptnr_idx].index;
            const index_type idx_b            = idx_buffer_map_[j];
            const real_type  pots_sum_b       = pots_sum_b_[idx_b];
            const real_type  pots_sum_b_coefa = pots_sum_b*inv_coefa_;

            if(std::max(pots_sum_a_coefb, pots_sum_b_coefa) <= 1.0){ continue; }

            const real_type pot_buff_ab = pot_buff_range[ptnr_idx];
            if(pots_sum_a_coefb < pots_sum_b_coefa)
            {
                pots_partial_sum_b_[idx_b] += pot_buff_ab;
            }
            else // pots_sum_a / coefb >= pots_sum_b / coefa
            {
                pots_partial_sum_a += pot_buff_ab;
            }
        }
    }

    // force calculation
    for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
    {
        const index_type i = participants_a[idx_a];

        const auto& partner               = partner_buffer_a_    [idx_a];
        const auto& potential_buff_range  = pot_buffer_range_a_  [idx_a];
        const auto& derivative_buff_range = deriv_buffer_range_a_[idx_a];

        const real_type        pots_sum_a         = pots_sum_a_        [idx_a];
        const coordinate_type& derivs_sum_a       = derivs_sum_a_      [idx_a];
        const real_type        pots_partial_sum_a = pots_partial_sum_a_[idx_a];
        const real_type        inv_pots_sum_a     = 1.0 / pots_sum_a;
        const real_type        coefb_pots_sum_a   = coefa_ * inv_pots_sum_a;
        const real_type        ep_coefb_pots_sum_a  = epsilon_ * coefb_pots_sum_a;
        const real_type        ep_coefb_pots_sum_a2 = ep_coefb_pots_sum_a * inv_pots_sum_a;

        for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
        {
            const index_type       j             = partner[ptnr_idx].index;
            const index_type       idx_b         = idx_buffer_map_[j];
            const real_type        pot_buff_ab   = potential_buff_range [ptnr_idx];
            const coordinate_type& deriv_buff_ab = derivative_buff_range[ptnr_idx];
            const real_type        pots_sum_b    = pots_sum_b_  [idx_b];
            const real_type        inv_pots_sum_b   = 1.0 / pots_sum_b;
            const real_type        coefa_pots_sum_b = coefa_ * inv_pots_sum_b;
            const coordinate_type& derivs_sum_b     = derivs_sum_b_[idx_b];

            if(1.0 <= coefb_pots_sum_a) // pots_sum_a / coefb <= 1.0 case
            {
                if(1.0 <= coefa_pots_sum_b) // pots_sum_b / coefa <= 1.0 case
                {
                    // force from distance part dirivation
                    const coordinate_type force = epsilon_ * deriv_buff_ab;
                    sys.force(i) -= force;
                    sys.force(j) += force;
                }
                else // pots_sum_a / coefb < pots_sum_b / coefa case
                {
                    const real_type ep_coefa_pots_sum_b  = epsilon_ * coefa_pots_sum_b;
                    const real_type ep_coefa_pots_sum_b2 = ep_coefa_pots_sum_b * inv_pots_sum_b;

                    // force from distance part derivation
                    sys.force(i) -= ep_coefa_pots_sum_b * deriv_buff_ab;
                    sys.force(j) += ep_coefa_pots_sum_b * deriv_buff_ab;

                    // force from valence part derivation
                    sys.force(j) += ep_coefa_pots_sum_b2 * derivs_sum_b * pot_buff_ab;
                }
            }
            else // 1.0 < pots_sum_a / coefb case
            {
                // force from depends on valence part
                if(coefb_pots_sum_a <= coefa_pots_sum_b) // pots_sum_b / coefa <= pots_sum_a / coefb
                {
                    // force from distance part derivation
                    sys.force(i) -= ep_coefb_pots_sum_a * deriv_buff_ab;
                    sys.force(j) += ep_coefb_pots_sum_a * deriv_buff_ab;

                    // force from valence part derivation
                    sys.force(i) += ep_coefb_pots_sum_a2 * derivs_sum_a * pot_buff_ab;
                }
                else // pots_sum_a / coefb < pots_sum_b / coefa
                {
                    const real_type ep_coefa_pots_sum_b  = epsilon_ * coefa_pots_sum_b;
                    const real_type ep_coefa_pots_sum_b2 = ep_coefa_pots_sum_b * inv_pots_sum_b;

                    // force from distance part derivation
                    sys.force(i) -= ep_coefa_pots_sum_b * deriv_buff_ab;
                    sys.force(j) += ep_coefa_pots_sum_b * deriv_buff_ab;
                    // force from valence part derivation
                    sys.force(j) += ep_coefa_pots_sum_b2 * derivs_sum_b * pot_buff_ab;
                }
                // force from depends on valence part
                sys.force(j) -= ep_coefb_pots_sum_a2 * deriv_buff_ab * pots_partial_sum_a;
            }
        }

        // force from depends on valence part
        for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
        {
            const index_type k              = partner[ptnr_idx].index;
            const index_type idx_d          = idx_buffer_map_[k];
            const real_type  pots_sum_d     = pots_sum_b_[idx_d];
            const real_type  inv_pots_sum_d = 1.0 / pots_sum_d;
            const real_type  coefa_pots_sum_d = coefa_ * inv_pots_sum_d;
            if(coefa_pots_sum_d < 1.0) // 1.0 < pots_sum_d / coefa
            {
                const real_type        ep_coefa_pots_sum_d  = epsilon_ * coefa_pots_sum_d;
                const real_type        ep_coefa_pots_sum_d2 = ep_coefa_pots_sum_d * inv_pots_sum_d;
                const coordinate_type& deriv_buff_ad        = derivative_buff_range[ptnr_idx];

                sys.force(i) += ep_coefa_pots_sum_d2 * deriv_buff_ad * pots_partial_sum_b_[idx_d];
            }
        }
    }

    // potential calculation
    real_type retval(0.0);
    for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a) 
    {
        const partner_range_type&   partner     = partner_buffer_a_  [idx_a];
        const potential_range_type& pots_buff_a = pot_buffer_range_a_[idx_a];
        for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
        {
            const index_type j            = partner[ptnr_idx].index;
            const index_type idx_b        = idx_buffer_map_[j];
            const real_type  pots_buff_ab = pots_buff_a[ptnr_idx];
            const real_type  max_normalized_contact_order =
                std::max(pots_sum_a_[idx_a]*inv_coefb_, pots_sum_b_[idx_b]*inv_coefa_);
            if(max_normalized_contact_order <= 1.0)
            {
                retval += pots_buff_ab;
            }
            else
            {
                retval += pots_buff_ab / max_normalized_contact_order;
            }
        }
    }
    retval *= -epsilon_;

    return retval;
}
} // mjolnir

#endif /* MJOLNIR_FORCEFIELD_STOICHIOMETRIC_GLOBAL_STOICHIOMETRIC_INTERACTION_HPP */
