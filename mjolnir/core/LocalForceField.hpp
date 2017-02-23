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

    template<std::size_t N>
    struct interaction_index_pair
    {
        typedef LocalInteractionBase<traitsT, N>       interaction_type;
        typedef std::unique_ptr<interaction_type>      interaction_ptr;
        typedef std::array<std::size_t, N>             index_list;

        interaction_index_pair(interaction_ptr&& interact, index_list&& idxs)
            : interaction(std::forward<interaction_ptr>(interact)),
              indices(std::forward<index_list>(idxs))
        {}
        interaction_ptr interaction;
        index_list      indices;
    };

    template<std::size_t N>
    using interaction_ptr = typename interaction_index_pair<N>::interaction_ptr;
    template<std::size_t N>
    using index_list = typename interaction_index_pair<N>::index_list;

  public:

    LocalForceField() = default;
    ~LocalForceField() = default;
    LocalForceField(LocalForceField const&) = delete;
    LocalForceField(LocalForceField&&)      = default;
    LocalForceField& operator=(LocalForceField const&) = delete;
    LocalForceField& operator=(LocalForceField&&)      = default;

    void      calc_force(particle_container_type& pcon);
    real_type calc_energy(const particle_container_type& pcon) const;

    void emplace_2body(interaction_ptr<2>&& i, index_list<2>&& ps);
    void emplace_3body(interaction_ptr<3>&& i, index_list<3>&& ps);
    void emplace_4body(interaction_ptr<4>&& i, index_list<4>&& ps);

    void reset_parameter(const std::string&, const real_type);

  private:

    std::vector<interaction_index_pair<2>> body2_;
    std::vector<interaction_index_pair<3>> body3_;
    std::vector<interaction_index_pair<4>> body4_;
};

template<typename traitsT>
inline void LocalForceField<traitsT>::emplace_2body(
        interaction_ptr<2>&& i, index_list<2>&& ps)
{
    body2_.emplace_back(std::forward<interaction_ptr<2>>(i),
                        std::forward<index_list<2>>(ps));
    return;
}

template<typename traitsT>
inline void LocalForceField<traitsT>::emplace_3body(
        interaction_ptr<3>&& i, index_list<3>&& ps)
{
    body3_.emplace_back(std::forward<interaction_ptr<3>>(i),
                        std::forward<index_list<3>>(ps));
    return;
}

template<typename traitsT>
inline void LocalForceField<traitsT>::emplace_4body(
        interaction_ptr<4>&& i, index_list<4>&& ps)
{
    body4_.emplace_back(std::forward<interaction_ptr<4>>(i),
                        std::forward<index_list<4>>(ps));
    return;
}

template<typename traitsT>
void LocalForceField<traitsT>::calc_force(particle_container_type& pcon)
{
    for(auto iter = body2_.cbegin(); iter != body2_.cend(); ++iter)
        iter->interaction->calc_force(
                pcon[iter->indices[0]], pcon[iter->indices[1]]);

    for(auto iter = body3_.cbegin(); iter != body3_.cend(); ++iter)
        iter->interaction->calc_force(
                pcon[iter->indices[0]], pcon[iter->indices[1]],
                pcon[iter->indices[2]]);

    for(auto iter = body4_.cbegin(); iter != body4_.cend(); ++iter)
        iter->interaction->calc_force(
                pcon[iter->indices[0]], pcon[iter->indices[1]],
                pcon[iter->indices[2]], pcon[iter->indices[3]]);
    return;
}

template<typename traitsT>
typename LocalForceField<traitsT>::real_type
LocalForceField<traitsT>::calc_energy(const particle_container_type& pcon) const
{
    real_type energy = 0.0;

    for(auto iter = body2_.cbegin(); iter != body2_.cend(); ++iter)
        energy += iter->interaction->calc_energy(
            pcon[iter->indices[0]], pcon[iter->indices[1]]);

    for(auto iter = body3_.cbegin(); iter != body3_.cend(); ++iter)
        energy += iter->interaction->calc_energy(
            pcon[iter->indices[0]], pcon[iter->indices[1]],
            pcon[iter->indices[2]]);

    for(auto iter = body4_.cbegin(); iter != body4_.cend(); ++iter)
        energy += iter->interaction->calc_energy(
            pcon[iter->indices[0]], pcon[iter->indices[1]],
            pcon[iter->indices[2]], pcon[iter->indices[3]]);

    return energy;
}

template<typename traitsT>
void LocalForceField<traitsT>::reset_parameter(
        const std::string& name, const real_type val)
{
    for(auto iter = body2_.begin(); iter != body2_.end(); ++iter)
        iter->interaction->reset_parameter(name, val);
    for(auto iter = body3_.begin(); iter != body3_.end(); ++iter)
        iter->interaction->reset_parameter(name, val);
    for(auto iter = body4_.begin(); iter != body4_.end(); ++iter)
        iter->interaction->reset_parameter(name, val);
    return;
}

} // mjolnir
#endif /* MJOLNIR_LOCAL_FORCE_FIELD */
