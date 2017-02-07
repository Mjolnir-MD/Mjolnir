#ifndef MJOLNIR_LOCAL_FORCE_FIELD
#define MJOLNIR_LOCAL_FORCE_FIELD
#include "BondLengthInteraction.hpp"
#include "BondAngleInteraction.hpp"
#include "DihedralAngleInteraction.hpp"
#include "LocalPotentialBase.hpp"
#include "ParticleContainer.hpp"
#include <utility>
#include <vector>
#include <array>
#include <memory>

namespace mjolnir
{

template<typename traitsT>
class LocalForceField
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::time_type       time_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef Particle<coordinate_type>             particle_type;
    typedef ParticleContainer<traits_type>        particle_container_type;
    typedef LocalPotentialBase<traits_type>       potential_base;
    typedef std::unique_ptr<potential_base>       potential_ptr;

    template<std::size_t N>
    struct interaction_potential_pair
    {
        typedef LocalInteractionBase<traitsT, N>         interaction_type;
        typedef std::unique_ptr<interaction_type>        interaction_ptr;
        typedef std::array<std::size_t, N>               index_list;
        typedef std::pair<index_list, potential_ptr>     potential_index_pair;
        typedef std::vector<potential_index_pair>        potential_array;

        interaction_potential_pair(interaction_ptr&& i, potential_array&& ps)
            : interaction(std::forward<interaction_ptr>(i)),
              potentials(std::forward<potential_array>(ps))
        {}
        interaction_ptr interaction;
        potential_array potentials;
    };

    template<std::size_t N>
    using interaction_ptr = typename interaction_potential_pair<N>::interaction_ptr;
    template<std::size_t N>
    using potential_array = typename interaction_potential_pair<N>::potential_array;

  public:

    LocalForceField() = default;
    ~LocalForceField() = default;
    LocalForceField(LocalForceField const&) = delete;
    LocalForceField(LocalForceField&&)      = default;
    LocalForceField& operator=(LocalForceField const&) = delete;
    LocalForceField& operator=(LocalForceField&&)      = default;

    void      calc_force(particle_container_type& pcon);
    real_type calc_energy(const particle_container_type& pcon) const;

    void emplace_2body(interaction_ptr<2>&& i, potential_array<2>&& ps);
    void emplace_3body(interaction_ptr<3>&& i, potential_array<3>&& ps);
    void emplace_4body(interaction_ptr<4>&& i, potential_array<4>&& ps);

  private:

    std::vector<interaction_potential_pair<2>> body2_;
    std::vector<interaction_potential_pair<3>> body3_;
    std::vector<interaction_potential_pair<4>> body4_;
};

template<typename traitsT>
inline void LocalForceField<traitsT>::emplace_2body(
        interaction_ptr<2>&& i, potential_array<2>&& ps)
{
    body2_.emplace_back(std::forward<interaction_ptr<2>>(i),
                        std::forward<potential_array<2>>(ps));
    return;
}

template<typename traitsT>
inline void LocalForceField<traitsT>::emplace_3body(
        interaction_ptr<3>&& i, potential_array<3>&& ps)
{
    body3_.emplace_back(std::forward<interaction_ptr<3>>(i),
                        std::forward<potential_array<3>>(ps));
    return;
}

template<typename traitsT>
inline void LocalForceField<traitsT>::emplace_4body(
        interaction_ptr<4>&& i, potential_array<4>&& ps)
{
    body4_.emplace_back(std::forward<interaction_ptr<4>>(i),
                        std::forward<potential_array<4>>(ps));
    return;
}

template<typename traitsT>
void LocalForceField<traitsT>::calc_force(particle_container_type& pcon)
{
    for(auto ff = body2_.cbegin(); ff != body2_.cend(); ++ff)
    {
        const auto& interaction = ff->interaction;
        const auto& potentials  = ff->potentials;
        for(auto iter = potentials.cbegin(); iter != potentials.cend(); ++iter)
        {
            interaction->calc_force(pcon[iter->first[0]], pcon[iter->first[1]],
                                    *(iter->second));
        }
    }

    for(auto ff = body3_.cbegin(); ff != body3_.cend(); ++ff)
    {
        const auto& interaction = ff->interaction;
        const auto& potentials  = ff->potentials;
        for(auto iter = potentials.cbegin(); iter != potentials.cend(); ++iter)
        {
            interaction->calc_force(pcon[iter->first[0]], pcon[iter->first[1]],
                                    pcon[iter->first[2]], *(iter->second));
        }
    }

    for(auto ff = body4_.cbegin(); ff != body4_.cend(); ++ff)
    {
        const auto& interaction = ff->interaction;
        const auto& potentials  = ff->potentials;
        for(auto iter = potentials.cbegin(); iter != potentials.cend(); ++iter)
            interaction->calc_force(pcon[iter->first[0]], pcon[iter->first[1]],
                                    pcon[iter->first[2]], pcon[iter->first[3]],
                                    *(iter->second));
    }
    return;
}

template<typename traitsT>
typename LocalForceField<traitsT>::real_type
LocalForceField<traitsT>::calc_energy(const particle_container_type& pcon) const
{
    real_type energy = 0.0;

    for(auto ff = body2_.cbegin(); ff != body2_.cend(); ++ff)
    {
        const auto& interaction = ff->interaction;
        const auto& potentials  = ff->potentials;
        for(auto iter = potentials.cbegin(); iter != potentials.cend(); ++iter)
            energy += interaction->calc_energy(
                pcon[iter->first[0]], pcon[iter->first[1]], *(iter->second));
    }

    for(auto ff = body3_.cbegin(); ff != body3_.cend(); ++ff)
    {
        const auto& interaction = ff->interaction;
        const auto& potentials  = ff->potentials;
        for(auto iter = potentials.cbegin(); iter != potentials.cend(); ++iter)
            energy += interaction->calc_energy(
                pcon[iter->first[0]], pcon[iter->first[1]], pcon[iter->first[2]],
                *(iter->second));
    }

    for(auto ff = body4_.cbegin(); ff != body4_.cend(); ++ff)
    {
        const auto& interaction = ff->interaction;
        const auto& potentials  = ff->potentials;
        for(auto iter = potentials.cbegin(); iter != potentials.cend(); ++iter)
            energy += interaction->calc_energy(
                pcon[iter->first[0]], pcon[iter->first[1]], pcon[iter->first[2]],
                pcon[iter->first[3]], *(iter->second));
    }
    return energy;
}

} // mjolnir
#endif /* MJOLNIR_LOCAL_FORCE_FIELD */
