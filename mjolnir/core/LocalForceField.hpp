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

  public:

    LocalForceField() = default;
    ~LocalForceField() = default;
    LocalForceField(LocalForceField const&) = delete;
    LocalForceField(LocalForceField&&)      = default;
    LocalForceField& operator=(LocalForceField const&) = delete;
    LocalForceField& operator=(LocalForceField&&)      = default;

    void initialize(const particle_container_type& pcon, const time_type dt);

    void emplace_bond(std::size_t i, std::size_t j, potential_ptr&& pot);
    void emplace_angle(std::size_t i, std::size_t j, std::size_t k,
                       potential_ptr&& pot);
    void emplace_dihedral(std::size_t i, std::size_t j, std::size_t k,
                          std::size_t l, potential_ptr&& pot);

    void      calc_force(particle_container_type& pcon);
    real_type calc_energy(const particle_container_type& pcon) const;

  private:

    BondLengthInteraction<traitsT>    bond;
    BondAngleInteraction<traitsT>     angle;
    DihedralAngleInteraction<traitsT> dihd;

    struct bond_potential_type
    {
        typedef LocalInteractionBase<traitsT, 2> interaction_type;
        typedef typename interaction_type::particle_ptrs particle_ptrs;
        bond_potential_type() = default;
        bond_potential_type(std::size_t i_, std::size_t j_, potential_ptr&& pot_)
           : indices{{i_, j_}}, pot(std::forward<potential_ptr>(pot_))
        {}
        ~bond_potential_type() = default;
        bond_potential_type(const bond_potential_type&) = delete;
        bond_potential_type(bond_potential_type&&)      = default;
        bond_potential_type& operator=(const bond_potential_type&) = delete;
        bond_potential_type& operator=(bond_potential_type&&)      = default;

        std::array<std::size_t, 2> indices;
        particle_ptrs              particles;
        potential_ptr              pot;
    };

    struct angle_potential_type
    {
        typedef LocalInteractionBase<traitsT, 3> interaction_type;
        typedef typename interaction_type::particle_ptrs particle_ptrs;
        angle_potential_type() = default;
        angle_potential_type(std::size_t i_, std::size_t j_, std::size_t k_,
                             potential_ptr&& pot_)
           : indices{{i_, j_, k_}}, pot(std::forward<potential_ptr>(pot_))
        {}
        ~angle_potential_type() = default;
        angle_potential_type(const angle_potential_type&) = delete;
        angle_potential_type(angle_potential_type&&)      = default;
        angle_potential_type& operator=(const angle_potential_type&) = delete;
        angle_potential_type& operator=(angle_potential_type&&)      = default;

        std::array<std::size_t, 3> indices;
        particle_ptrs              particles;
        potential_ptr              pot;
    };

    struct dihd_potential_type
    {
        typedef LocalInteractionBase<traitsT, 4> interaction_type;
        typedef typename interaction_type::particle_ptrs particle_ptrs;
        dihd_potential_type() = default;
        dihd_potential_type(std::size_t i_, std::size_t j_, std::size_t k_,
                std::size_t l_, potential_ptr&& pot_)
           : indices{{i_, j_, k_, l_}}, pot(std::forward<potential_ptr>(pot_))
        {}
        ~dihd_potential_type() = default;
        dihd_potential_type(const dihd_potential_type&) = delete;
        dihd_potential_type(dihd_potential_type&&)      = default;
        dihd_potential_type& operator=(const dihd_potential_type&) = delete;
        dihd_potential_type& operator=(dihd_potential_type&&)      = default;

        std::array<std::size_t, 4> indices;
        particle_ptrs              particles;
        potential_ptr              pot;
    };

    std::vector<bond_potential_type>  bond_potentials;
    std::vector<angle_potential_type> angle_potentials;
    std::vector<dihd_potential_type>  dihd_potentials;
};

template<typename traitsT>
void LocalForceField<traitsT>::emplace_bond(
        std::size_t i, std::size_t j, potential_ptr&& pot)
{
    bond_potentials.emplace_back(
            i, j, std::forward<potential_ptr>(pot));
    return ;
}

template<typename traitsT>
void LocalForceField<traitsT>::emplace_angle(
        std::size_t i, std::size_t j, std::size_t k, potential_ptr&& pot)
{
    angle_potentials.emplace_back(i, j, k,
            std::forward<potential_ptr>(pot));
    return ;
}

template<typename traitsT>
void LocalForceField<traitsT>::emplace_dihedral(
        std::size_t i, std::size_t j, std::size_t k, std::size_t l,
        potential_ptr&& pot)
{
    dihd_potentials.emplace_back(i, j, k, l,
            std::forward<potential_ptr>(pot));
    return ;
}

template<typename traitsT>
void LocalForceField<traitsT>::initialize(
        const particle_container_type& pcon, const time_type dt)
{
    for(auto iter = bond_potentials.begin(); iter != bond_potentials.end(); ++iter)
    {
        iter->particles[0] = const_cast<particle_type*>(&(pcon.at(iter->indices[0])));
        iter->particles[1] = const_cast<particle_type*>(&(pcon.at(iter->indices[1])));
    }
    for(auto iter = angle_potentials.begin(); iter != angle_potentials.end(); ++iter)
    {
        iter->particles[0] = const_cast<particle_type*>(&(pcon.at(iter->indices[0])));
        iter->particles[1] = const_cast<particle_type*>(&(pcon.at(iter->indices[1])));
        iter->particles[2] = const_cast<particle_type*>(&(pcon.at(iter->indices[2])));
    }
    for(auto iter = dihd_potentials.begin(); iter != dihd_potentials.end(); ++iter)
    {
        iter->particles[0] = const_cast<particle_type*>(&(pcon.at(iter->indices[0])));
        iter->particles[1] = const_cast<particle_type*>(&(pcon.at(iter->indices[1])));
        iter->particles[2] = const_cast<particle_type*>(&(pcon.at(iter->indices[2])));
        iter->particles[3] = const_cast<particle_type*>(&(pcon.at(iter->indices[3])));
    }
    return;
}

template<typename traitsT>
void LocalForceField<traitsT>::calc_force(particle_container_type& pcon)
{
    for(auto iter = bond_potentials.cbegin();
            iter != bond_potentials.cend(); ++iter)
        bond.calc_force(iter->particles, *(iter->pot));

    for(auto iter = angle_potentials.cbegin();
            iter != angle_potentials.cend(); ++iter)
        angle.calc_force(iter->particles, *(iter->pot));

    for(auto iter = dihd_potentials.cbegin();
            iter != dihd_potentials.cend(); ++iter)
        dihd.calc_force(iter->particles, *(iter->pot));

    return;
}

template<typename traitsT>
typename LocalForceField<traitsT>::real_type
LocalForceField<traitsT>::calc_energy(const particle_container_type& pcon) const
{
    real_type energy = 0.0;

    for(auto iter = bond_potentials.cbegin();
            iter != bond_potentials.cend(); ++iter)
        energy += bond.calc_energy(iter->particles, *(iter->pot));

    for(auto iter = angle_potentials.cbegin();
            iter != angle_potentials.cend(); ++iter)
        energy += angle.calc_energy(iter->particles, *(iter->pot));

    for(auto iter = dihd_potentials.cbegin();
            iter != dihd_potentials.cend(); ++iter)
        energy += dihd.calc_energy(iter->particles, *(iter->pot));
    return energy;
}

} // mjolnir
#endif /* MJOLNIR_LOCAL_FORCE_FIELD */
