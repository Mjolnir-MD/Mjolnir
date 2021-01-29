#ifndef MJOLNIR_OMP_GLOBAL_STOICHIOMETRIC_INTERACTION_HPP
#define MJOLNIR_OMP_GLOBAL_STOICHIOMETRIC_INTERACTION_HPP
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/forcefield/stoichiometric/GlobalStoichiometricInteraction.hpp>

namespace mjolnir
{

// The interaction which conserve stoichiometry.

template<typename realT, template<typename, typename> class boundaryT>
class GlobalStoichiometricInteraction<OpenMPSimulatorTraits<realT, boundaryT>>
    final : public GlobalInteractionBase<OpenMPSimulatorTraits<realT, boundaryT>>
{
  public:

    using traits_type            = OpenMPSimulatorTraits<realT, boundaryT>;
    using base_type              = GlobalInteractionBase<traits_type>;
    using index_type             = std::size_t;
    using real_type              = typename base_type::real_type;
    using coordinate_type        = typename base_type::coordinate_type;
    using system_type            = typename base_type::system_type;
    using topology_type          = typename base_type::topology_type;
    using potential_type         = GlobalStoichiometricInteractionPotential<traits_type>;
    using partition_type         = SpatialPartition<traits_type, potential_type>;
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
        num_threads_ = omp_get_max_threads();

        const std::size_t participants_a_num = potential_.participants_a_num();
        const std::size_t participants_b_num = potential_.participants_b_num();
        const std::size_t participants_ab    = participants_a_num * participants_b_num;
        pot_derivs_buff_.resize(participants_ab);
        potentials_buff_.resize(participants_ab);
        partner_buffer_ranges_.resize(participants_a_num);

        pot_sum_a_      .resize(participants_a_num);
        pot_sum_b_      .resize(participants_b_num);
        pot_deriv_sum_a_.resize(participants_a_num);
        pot_deriv_sum_b_.resize(participants_b_num);

        pot_sum_thread_b_      .resize(num_threads_);
        pot_deriv_sum_thread_b_.resize(num_threads_);
        for(std::size_t thread_idx=0; thread_idx<num_threads_; ++thread_idx)
        {
            pot_sum_thread_b_[thread_idx]      .resize(participants_b_num);
            pot_deriv_sum_thread_b_[thread_idx].resize(participants_b_num);
        }
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
        idx_buffer_map.resize(sys.size());
        const auto following_participants = potential_.following_participants();
        for(std::size_t idx_b=0; idx_b<following_participants.size(); ++idx_b)
        {
            idx_buffer_map[following_participants[idx_b]] = idx_b;
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
        const auto        leading_participants = potential_.leading_participants();
        const std::size_t participants_a_num   = potential_.participants_a_num();
        auto pot_buff_iter   = potentials_buff_.begin();
        auto deriv_buff_iter = pot_derivs_buff_.begin();

        for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
        {
            const index_type i = leading_participants[idx_a];

            const partner_range_type& partner    = partition_.partners(i);
            const index_type          range_size = partner.size();
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
    std::size_t    num_threads_;

    // -----------------------------------------------------------------------
    // Variable for mapping system index to buffering index
    std::vector<std::size_t> idx_buffer_map;

    // Variables for buffering intermediate value.
    // These value can change in calc_force and calc_energy function.
    mutable potential_buffer_type  potentials_buff_;
    mutable derivative_buffer_type pot_derivs_buff_;
    mutable potential_buffer_type  pot_sum_a_;
    mutable potential_buffer_type  pot_sum_b_;

    // Variables to specify the range of specific a partners in buffer.
    mutable partner_buffer_ranges_type partner_buffer_ranges_;

    // To avoid data race during summation about b
    mutable std::vector<potential_buffer_type>  pot_sum_thread_b_;
    mutable std::vector<derivative_buffer_type> pot_deriv_sum_thread_b_;

    // sum of derivation of potential function for specific first particle.
    mutable derivative_buffer_type pot_deriv_sum_a_;

    // sum of derivation of potential function for specific second particle.
    mutable derivative_buffer_type pot_deriv_sum_b_;

  public:

    void calc_force(system_type& sys) const noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
        MJOLNIR_LOG_FUNCTION_DEBUG();
    
        // initialization of each buffering container.
        std::fill(pot_sum_a_.begin(), pot_sum_a_.end(), 0.0);
        std::fill(pot_sum_b_.begin(), pot_sum_b_.end(), 0.0);
        std::fill(pot_deriv_sum_a_.begin(), pot_deriv_sum_a_.end(),
                  math::make_coordinate<coordinate_type>(0.0, 0.0, 0.0));
        std::fill(pot_deriv_sum_b_.begin(), pot_deriv_sum_b_.end(),
                  math::make_coordinate<coordinate_type>(0.0, 0.0, 0.0));
        for(std::size_t thread_idx=0; thread_idx<num_threads_; ++thread_idx)
        {
            potential_buffer_type&  pot_sum_thread       = pot_sum_thread_b_      [thread_idx];
            derivative_buffer_type& pot_deriv_sum_thread = pot_deriv_sum_thread_b_[thread_idx];
            std::fill(pot_sum_thread.begin(),       pot_sum_thread.end(), 0.0);
            std::fill(pot_deriv_sum_thread.begin(), pot_deriv_sum_thread.end(),
                      math::make_coordinate<coordinate_type>(0.0, 0.0, 0.0));
        }

        // make index pair list to specify the range of buffer for specific a and
        // pre calculation for pot_sum and pot_deriv_sum for each particle.
        const auto        leading_participants = potential_.leading_participants();
        const std::size_t participants_a_num   = potential_.participants_a_num();
        const std::size_t participants_b_num   = potential_.participants_b_num();

#pragma omp parallel for
        for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
        {
            const index_type i = leading_participants[idx_a];

            partner_buffer_tuple_type& partner_buffer_range = partner_buffer_ranges_[idx_a];
            const auto& partner     = std::get<0>(partner_buffer_range);
            const auto& pot_range   = std::get<1>(partner_buffer_range);
            const auto& deriv_range = std::get<2>(partner_buffer_range);

            real_type&       pot_sum_a       = pot_sum_a_      [idx_a];
            coordinate_type& pot_deriv_sum_a = pot_deriv_sum_a_[idx_a];

            const std::size_t thread_id = omp_get_thread_num();
            potential_buffer_type&  pot_sum_thread       = pot_sum_thread_b_      [thread_id];
            derivative_buffer_type& pot_deriv_sum_thread = pot_deriv_sum_thread_b_[thread_id];

            for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
            {
                const index_type j     = partner[ptnr_idx].index;
                const coordinate_type rij =
                    sys.adjust_direction(sys.position(i), sys.position(j));
                const real_type       l2    = math::length_sq(rij); // |rij|^2
                const real_type       rl    = math::rsqrt(l2);      // 1 / |rij|
                const real_type       l     = l2 * rl;
                const coordinate_type deriv = potential_.derivative(l) * rl * rij;
                const real_type       pot   = potential_.potential(l);

                pot_range        [ptnr_idx] =  pot;
                deriv_range      [ptnr_idx] =  deriv;
                pot_sum_a                   += pot;
                pot_deriv_sum_a             += deriv;

                const index_type     idx_b  =  idx_buffer_map[j];
                pot_sum_thread      [idx_b] += pot;
                pot_deriv_sum_thread[idx_b] -= deriv;
            }
        }

        for(std::size_t thread_idx=0; thread_idx<num_threads_; ++thread_idx)
        {
            const potential_buffer_type&  pot_sum_thread       = pot_sum_thread_b_      [thread_idx];
            const derivative_buffer_type& pot_deriv_sum_thread = pot_deriv_sum_thread_b_[thread_idx];
            for(std::size_t idx_b=0; idx_b<participants_b_num; ++idx_b)
            {
                pot_sum_b_[idx_b]       += pot_sum_thread      [idx_b];
                pot_deriv_sum_b_[idx_b] += pot_deriv_sum_thread[idx_b];
            }
        }

        // force calculation
#pragma omp parallel for
        for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
        {
            const index_type i = leading_participants[idx_a];

            const auto&     partner_buffer_tuple  = partner_buffer_ranges_[idx_a];
            const auto&     partner               = std::get<0>(partner_buffer_tuple);
            const auto&     potential_buff_range  = std::get<1>(partner_buffer_tuple);
            const auto&     derivative_buff_range = std::get<2>(partner_buffer_tuple);
            const real_type pot_sum_a             = pot_sum_a_[idx_a];

            const std::size_t thread_id = omp_get_thread_num();
    
            if(pot_sum_a <= 1.0)
            {
                for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
                {
                    const index_type       j         = partner[ptnr_idx].index;
                    const index_type       idx_b     = idx_buffer_map[j];
                    const real_type        pot_sum_b = pot_sum_b_[idx_b];
                    const coordinate_type& pot_derivs_buff_ab = derivative_buff_range[ptnr_idx];
    
                    if(pot_sum_b <= 1.0)
                    {
                        const coordinate_type force = epsilon_ * pot_derivs_buff_ab;
                        sys.force_thread(thread_id, i) -= force;
                        sys.force_thread(thread_id, j) += force;
                    } 
                    else
                    {
                        const real_type pots_buff_ab           = potential_buff_range[ptnr_idx];
                        const real_type pot_sum_b              = pot_sum_b_             [idx_b];
                        const real_type inv_pot_sum_b          = 1.0 / pot_sum_b;
                        const real_type ep_pot_sum_b           = epsilon_ * inv_pot_sum_b;
                        const real_type pots_buff_ab_pot_sum_b = pots_buff_ab * inv_pot_sum_b;
                        sys.force_thread(thread_id, i) +=
                            ep_pot_sum_b *
                            (pots_buff_ab_pot_sum_b - 1.0) * pot_derivs_buff_ab;
                        sys.force_thread(thread_id, j) -=
                            ep_pot_sum_b *
                            (pots_buff_ab_pot_sum_b * pot_deriv_sum_b_[idx_b] - pot_derivs_buff_ab);
                    }
                }
            }
            else // 1.0 < pot_sum_a case
            {
                const real_type inv_pot_sum_a  = 1.0 / pot_sum_a;
                const real_type ep_pot_sum_a   = epsilon_ * inv_pot_sum_a;
                const coordinate_type& pot_derivs_sum_a = pot_deriv_sum_a_[idx_a];
                const coordinate_type  pot_derivs_sum_a_pot_sum_a = inv_pot_sum_a * pot_derivs_sum_a;
    
                for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
                {
                    const index_type j     = partner[ptnr_idx].index;
                    const index_type idx_b = idx_buffer_map[j];
                    const real_type        pots_buff_ab       = potential_buff_range [ptnr_idx];
                    const coordinate_type& pot_derivs_buff_ab = derivative_buff_range[ptnr_idx];
                    const real_type        pot_sum_b          = pot_sum_b_           [idx_b];
    
                    if(pot_sum_b < pot_sum_a)
                    {
                        sys.force_thread(thread_id, i) +=
                            ep_pot_sum_a *
                            (pots_buff_ab * pot_derivs_sum_a_pot_sum_a - pot_derivs_buff_ab);
                        sys.force_thread(thread_id, j) -=
                            ep_pot_sum_a *
                            (pots_buff_ab * inv_pot_sum_a - 1.0) * pot_derivs_buff_ab;
                    }
                    else // pot_sum_a <= pot_sum_b
                    {
                        const real_type inv_pot_sum_b  = 1.0 / pot_sum_b;
                        const real_type ep_pot_sum_b   = epsilon_ * inv_pot_sum_b;
                        const real_type pots_buff_ab_pot_sum_b = pots_buff_ab * inv_pot_sum_b;
                        sys.force_thread(thread_id, i) +=
                            ep_pot_sum_b *
                            (pots_buff_ab_pot_sum_b - 1.0) * pot_derivs_buff_ab;
                        sys.force_thread(thread_id, j) -=
                            ep_pot_sum_b *
                            (pots_buff_ab_pot_sum_b * pot_deriv_sum_b_[idx_b] - pot_derivs_buff_ab);
                    }
                }
            }
        }
        return;
    }
    
    real_type calc_energy(const system_type& sys) const noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
        MJOLNIR_LOG_FUNCTION_DEBUG();

        // initialization of each buffering container.
        std::fill(pot_sum_a_.begin(), pot_sum_a_.end(), 0.0);
        std::fill(pot_sum_b_.begin(), pot_sum_b_.end(), 0.0);
        for(std::size_t thread_idx=0; thread_idx<num_threads_; ++thread_idx)
        {
            potential_buffer_type& pot_sum_thread = pot_sum_thread_b_[thread_idx];
            std::fill(pot_sum_thread.begin(), pot_sum_thread.end(), 0.0);
        }

        // make index pair list to specify the range of buffer for specific a and
        // pre calculation for pot_sum and pot_deriv_sum for each particle.
        const auto        leading_participants = potential_.leading_participants();
        const std::size_t participants_a_num   = potential_.participants_a_num();
        const std::size_t participants_b_num   = potential_.participants_b_num();

#pragma omp parallel for
        for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
        {
            const index_type i = leading_participants[idx_a];

            partner_buffer_tuple_type& partner_buffer_range = partner_buffer_ranges_[idx_a];
            const auto& partner   = std::get<0>(partner_buffer_range);
            const auto& pot_range = std::get<1>(partner_buffer_range);

            real_type& pot_sum_a = pot_sum_a_[idx_a];

            const std::size_t thread_id = omp_get_thread_num();
            potential_buffer_type&  pot_sum_thread = pot_sum_thread_b_[thread_id];

            for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
            {
                const index_type j     = partner[ptnr_idx].index;
                const index_type idx_b = idx_buffer_map[j];
                const coordinate_type rij =
                    sys.adjust_direction(sys.position(i), sys.position(j));
                const real_type l2  = math::length_sq(rij); // |rij|^2
                const real_type rl  = math::rsqrt(l2);      // 1 / |rij|
                const real_type l   = l2 * rl;              // |rij|
                const real_type pot = potential_.potential(l);

                pot_range  [ptnr_idx] =  pot;
                pot_sum_a             += pot;
                pot_sum_thread[idx_b] += pot;
            }
        }

        for(std::size_t thread_idx=0; thread_idx<num_threads_; ++thread_idx)
        {
            const potential_buffer_type& pot_sum_thread = pot_sum_thread_b_[thread_idx];
            for(std::size_t idx_b=0; idx_b<participants_b_num; ++idx_b)
            {
                pot_sum_b_[idx_b] += pot_sum_thread[idx_b];
            }
        }

        // potential calculation
        real_type retval(0.0);
#pragma omp parallel for reduction(+:retval)
        for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a) 
        {
            const partner_buffer_tuple_type&   partner_buffer_range = partner_buffer_ranges_[idx_a];
            const partner_range_type&          partner              = std::get<0>(partner_buffer_range);
            const potential_buffer_range_type& pots_buff_a          = std::get<1>(partner_buffer_range);
            for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
            {
                const index_type j            = partner[ptnr_idx].index;
                const index_type idx_b        = idx_buffer_map[j];
                const real_type  pots_buff_ab = pots_buff_a[ptnr_idx];
                const real_type  max_pot_sum  = std::max(pot_sum_a_[idx_a], pot_sum_b_[idx_b]);
                if(max_pot_sum <= 1.0)
                {
                    retval += pots_buff_ab;
                }
                else
                {
                    retval += pots_buff_ab / max_pot_sum;
                }
            }
        }
        retval *= -epsilon_;
    
        return retval;
    }
    
    real_type calc_force_and_energy(system_type& sys) const noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
        MJOLNIR_LOG_FUNCTION_DEBUG();
    
        // initialization of each buffering container.
        std::fill(pot_sum_a_.begin(), pot_sum_a_.end(), 0.0);
        std::fill(pot_sum_b_.begin(), pot_sum_b_.end(), 0.0);
        std::fill(pot_deriv_sum_a_.begin(), pot_deriv_sum_a_.end(),
                  math::make_coordinate<coordinate_type>(0.0, 0.0, 0.0));
        std::fill(pot_deriv_sum_b_.begin(), pot_deriv_sum_b_.end(),
                  math::make_coordinate<coordinate_type>(0.0, 0.0, 0.0));
        for(std::size_t thread_idx=0; thread_idx<num_threads_; ++thread_idx)
        {
            potential_buffer_type&  pot_sum_thread       = pot_sum_thread_b_      [thread_idx];
            derivative_buffer_type& pot_deriv_sum_thread = pot_deriv_sum_thread_b_[thread_idx];
            std::fill(pot_sum_thread.begin(),       pot_sum_thread.end(), 0.0);
            std::fill(pot_deriv_sum_thread.begin(), pot_deriv_sum_thread.end(),
                      math::make_coordinate<coordinate_type>(0.0, 0.0, 0.0));
        }

        // make index pair list to specify the range of buffer for specific a and
        // pre calculation for pot_sum and pot_deriv_sumb for each particle.
        const auto        leading_participants = potential_.leading_participants();
        const std::size_t participants_a_num   = potential_.participants_a_num();
        const std::size_t participants_b_num   = potential_.participants_b_num();

#pragma omp parallel for
        for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
        {
            const index_type i = leading_participants[idx_a];

            partner_buffer_tuple_type& partner_buffer_range = partner_buffer_ranges_[idx_a];
            const auto& partner     = std::get<0>(partner_buffer_range);
            const auto& pot_range   = std::get<1>(partner_buffer_range);
            const auto& deriv_range = std::get<2>(partner_buffer_range);

            real_type&       pot_sum_a       = pot_sum_a_      [idx_a];
            coordinate_type& pot_deriv_sum_a = pot_deriv_sum_a_[idx_a];

            const std::size_t thread_id = omp_get_thread_num();
            potential_buffer_type&  pot_sum_thread       = pot_sum_thread_b_      [thread_id];
            derivative_buffer_type& pot_deriv_sum_thread = pot_deriv_sum_thread_b_[thread_id];
    
            for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
            {
                const index_type j     = partner[ptnr_idx].index;
                const index_type idx_b = idx_buffer_map[j];
                const coordinate_type rij =
                   sys.adjust_direction(sys.position(i), sys.position(j));
                const real_type       l2    = math::length_sq(rij); // |rij|^2
                const real_type       rl    = math::rsqrt(l2);      // 1 / |rij|
                const real_type       l     = l2 * rl;
                const coordinate_type deriv = potential_.derivative(l) * rl * rij;
                const real_type       pot   = potential_.potential(l);

                pot_range        [ptnr_idx] =  pot;
                deriv_range      [ptnr_idx] =  deriv;
                pot_sum_a                   += pot;
                pot_deriv_sum_a             += deriv;
                pot_sum_thread      [idx_b] += pot;
                pot_deriv_sum_thread[idx_b] -= deriv;
            }
        }
    
        for(std::size_t thread_idx=0; thread_idx<num_threads_; ++thread_idx)
        {
            const potential_buffer_type&  pot_sum_thread       = pot_sum_thread_b_      [thread_idx];
            const derivative_buffer_type& pot_deriv_sum_thread = pot_deriv_sum_thread_b_[thread_idx];
            for(std::size_t idx_b=0; idx_b<participants_b_num; ++idx_b)
            {
                pot_sum_b_[idx_b]       += pot_sum_thread      [idx_b];
                pot_deriv_sum_b_[idx_b] += pot_deriv_sum_thread[idx_b];
            }
        }
    
        // force calculation
#pragma omp parallel for
        for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
        {
            const index_type i = leading_participants[idx_a];

            const auto&     partner_buffer_tuple  = partner_buffer_ranges_[idx_a];
            const auto&     partner               = std::get<0>(partner_buffer_tuple);
            const auto&     potential_buff_range  = std::get<1>(partner_buffer_tuple);
            const auto&     derivative_buff_range = std::get<2>(partner_buffer_tuple);
            const real_type pot_sum_a             = pot_sum_a_[idx_a];

            const std::size_t thread_id = omp_get_thread_num();

            if(pot_sum_a <= 1.0)
            {
                for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
                {
                    const index_type       j         = partner[ptnr_idx].index;
                    const index_type       idx_b     = idx_buffer_map[j];
                    const real_type        pot_sum_b = pot_sum_b_[idx_b];
                    const coordinate_type& pot_derivs_buff_ab = derivative_buff_range[ptnr_idx];
    
                    if(pot_sum_b <= 1.0)
                    {
                        const coordinate_type force = epsilon_ * pot_derivs_buff_ab;
                        sys.force_thread(thread_id, i) -= force;
                        sys.force_thread(thread_id, j) += force;
                    }
                    else
                    {
                        const real_type pots_buff_ab           = potential_buff_range[ptnr_idx];
                        const real_type pot_sum_b              = pot_sum_b_             [idx_b];
                        const real_type inv_pot_sum_b          = 1.0 / pot_sum_b;
                        const real_type ep_pot_sum_b           = epsilon_ * inv_pot_sum_b;
                        const real_type pots_buff_ab_pot_sum_b = pots_buff_ab * inv_pot_sum_b;
                        sys.force_thread(thread_id, i) +=
                            ep_pot_sum_b *
                            (pots_buff_ab_pot_sum_b - 1.0) * pot_derivs_buff_ab;
                        sys.force_thread(thread_id, j) -=
                            ep_pot_sum_b *
                            (pots_buff_ab_pot_sum_b * pot_deriv_sum_b_[idx_b] - pot_derivs_buff_ab);
                    }
                }
            }
            else // 1.0 < pot_sum_a case
            {
                const real_type inv_pot_sum_a = 1.0 / pot_sum_a;
                const real_type ep_pot_sum_a  = epsilon_ * inv_pot_sum_a;
                const coordinate_type& pot_derivs_sum_a = pot_deriv_sum_a_[idx_a];
                const coordinate_type  pot_derivs_sum_a_pot_sum_a = inv_pot_sum_a * pot_derivs_sum_a;
    
                for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
                {
                    const index_type j     = partner[ptnr_idx].index;
                    const index_type idx_b = idx_buffer_map[j];
                    const real_type        pots_buff_ab       = potential_buff_range [ptnr_idx];
                    const coordinate_type& pot_derivs_buff_ab = derivative_buff_range[ptnr_idx];
                    const real_type        pot_sum_b          = pot_sum_b_           [idx_b];
    
                    if(pot_sum_b < pot_sum_a)
                    {
                        sys.force_thread(thread_id, i) +=
                            ep_pot_sum_a *
                            (pots_buff_ab * pot_derivs_sum_a_pot_sum_a - pot_derivs_buff_ab);
                        sys.force_thread(thread_id, j) -=
                            ep_pot_sum_a *
                            (pots_buff_ab * inv_pot_sum_a - 1.0) * pot_derivs_buff_ab;
                    }
                    else // pot_sum_a <= pot_sum_b
                    {
                        const real_type inv_pot_sum_b = 1.0 / pot_sum_b;
                        const real_type ep_pot_sum_b  = epsilon_ * inv_pot_sum_b;
                        const real_type pots_buff_ab_pot_sum_b = pots_buff_ab * inv_pot_sum_b;
                        sys.force_thread(thread_id, i) +=
                            ep_pot_sum_b *
                            (pots_buff_ab_pot_sum_b - 1.0) * pot_derivs_buff_ab;
                        sys.force_thread(thread_id, j) -=
                            ep_pot_sum_b *
                            (pots_buff_ab_pot_sum_b * pot_deriv_sum_b_[idx_b] - pot_derivs_buff_ab);
                    }
                }
            }
        }
    
        // potential calculation
        real_type retval(0.0);
#pragma omp parallel for reduction(+:retval)
        for(std::size_t idx_a=0; idx_a<participants_a_num; ++idx_a)
        {
            const partner_buffer_tuple_type&   partner_buffer_range = partner_buffer_ranges_[idx_a];
            const partner_range_type&          partner              = std::get<0>(partner_buffer_range);
            const potential_buffer_range_type& pots_buff_a          = std::get<1>(partner_buffer_range);
            for(std::size_t ptnr_idx=0; ptnr_idx<partner.size(); ++ptnr_idx)
            {
                const index_type j            = partner[ptnr_idx].index;
                const index_type idx_b        = idx_buffer_map[j];
                const real_type  pots_buff_ab = pots_buff_a[ptnr_idx];
                const real_type  max_pot_sum  = std::max(pot_sum_a_[idx_a], pot_sum_b_[idx_b]);
                if(max_pot_sum <= 1.0)
                {
                    retval += pots_buff_ab;
                }
                else
                {
                    retval += pots_buff_ab / max_pot_sum;
                }
            }
        }
        retval *= -epsilon_;
    
        return retval;
    }
};
} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize major use-cases
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class GlobalStoichiometricInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>       >;
extern template class GlobalStoichiometricInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class GlobalStoichiometricInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class GlobalStoichiometricInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_OMP_GLOBAL_STOICHIOMETRIC_INTERACTION_HPP */
