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
    struct local_forcefield_type
    {
        typedef LocalInteractionBase<traitsT, N>         interaction_type;
        typedef std::unique_ptr<interaction_type>        interaction_ptr;
        typedef typename interaction_type::particle_ptrs particle_ptrs;

        interaction_ptr interaction;
        std::vector<std::pair<particle_ptrs, potential_ptr>> potentials;
    };

  public:

    LocalForceField() = default;
    ~LocalForceField() = default;
    LocalForceField(LocalForceField const&) = delete;
    LocalForceField(LocalForceField&&)      = default;
    LocalForceField& operator=(LocalForceField const&) = delete;
    LocalForceField& operator=(LocalForceField&&)      = default;

//     void initialize(const particle_container_type& pcon, const time_type dt);

    void      calc_force(particle_container_type& pcon);
    real_type calc_energy(const particle_container_type& pcon) const;

//     void emplace_bond(std::size_t i, std::size_t j, potential_ptr&& pot);
//     void emplace_angle(std::size_t i, std::size_t j, std::size_t k,
//                        potential_ptr&& pot);
//     void emplace_dihedral(std::size_t i, std::size_t j, std::size_t k,
//                           std::size_t l, potential_ptr&& pot);

    local_forcefield_type<2>& body2_at(std::size_t i){return body2.at(i);}
    local_forcefield_type<3>& body3_at(std::size_t i){return body3.at(i);}
    local_forcefield_type<4>& body4_at(std::size_t i){return body4.at(i);}

  private:

    std::vector<local_forcefield_type<2>> body2;
    std::vector<local_forcefield_type<3>> body3;
    std::vector<local_forcefield_type<4>> body4;
};

// template<typename traitsT>
// void LocalForceField<traitsT>::emplace_bond(
//         std::size_t i, std::size_t j, potential_ptr&& pot)
// {
//     bond_potentials.emplace_back(
//             i, j, std::forward<potential_ptr>(pot));
//     return ;
// }
//
// template<typename traitsT>
// void LocalForceField<traitsT>::emplace_angle(
//         std::size_t i, std::size_t j, std::size_t k, potential_ptr&& pot)
// {
//     angle_potentials.emplace_back(i, j, k,
//             std::forward<potential_ptr>(pot));
//     return ;
// }
//
// template<typename traitsT>
// void LocalForceField<traitsT>::emplace_dihedral(
//         std::size_t i, std::size_t j, std::size_t k, std::size_t l,
//         potential_ptr&& pot)
// {
//     dihd_potentials.emplace_back(i, j, k, l,
//             std::forward<potential_ptr>(pot));
//     return ;
// }

template<typename traitsT>
void LocalForceField<traitsT>::calc_force(particle_container_type& pcon)
{
    for(auto ff = body2.cbegin(); ff != body2.cend(); ++ff)
    {
        for(auto iter = ff->potentials.cbegin();
                iter != ff->potentials.cend(); ++iter)
            ff->interaction.calc_force(iter->first, *(iter->second));
    }

    for(auto ff = body3.cbegin(); ff != body3.cend(); ++ff)
    {
        for(auto iter = ff->potentials.cbegin();
                iter != ff->potentials.cend(); ++iter)
            ff->interaction.calc_force(iter->first, *(iter->second));
    }

    for(auto ff = body4.cbegin(); ff != body4.cend(); ++ff)
    {
        for(auto iter = ff->potentials.cbegin();
                iter != ff->potentials.cend(); ++iter)
            ff->interaction.calc_force(iter->first, *(iter->second));
    }
    return;
}

template<typename traitsT>
typename LocalForceField<traitsT>::real_type
LocalForceField<traitsT>::calc_energy(const particle_container_type& pcon) const
{
    real_type energy = 0.0;

    for(auto ff = body2.cbegin(); ff != body2.cend(); ++ff)
    {
        for(auto iter = ff->potentials.cbegin();
                iter != ff->potentials.cend(); ++iter)
            energy += ff->interaction.calc_energy(iter->first, *(iter->second));
    }

    for(auto ff = body3.cbegin(); ff != body3.cend(); ++ff)
    {
        for(auto iter = ff->potentials.cbegin();
                iter != ff->potentials.cend(); ++iter)
            energy += ff->interaction.calc_energy(iter->first, *(iter->second));
    }

    for(auto ff = body4.cbegin(); ff != body4.cend(); ++ff)
    {
        for(auto iter = ff->potentials.cbegin();
                iter != ff->potentials.cend(); ++iter)
            energy += ff->interaction.calc_energy(iter->first, *(iter->second));
    }
    return energy;
}

} // mjolnir
#endif /* MJOLNIR_LOCAL_FORCE_FIELD */
