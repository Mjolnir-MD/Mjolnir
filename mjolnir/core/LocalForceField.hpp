#ifndef MJOLNIR_LOCAL_FORCE_FIELD
#define MJOLNIR_LOCAL_FORCE_FIELD
#include "BondLengthInteraction.hpp"
#include "BondAngleInteraction.hpp"
#include "DihedralAngleInteraction.hpp"
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
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef Particle<coordinate_type> particle_type;
    typedef ParticleContainer<traits_type> particle_container_type;
    typedef LocalPotentialBase<traits_type> potential_base;
    typedef std::unique_ptr<potential_base> potential_ptr;
    typedef std::array<std::size_t, 2> pid2_type;
    typedef std::array<std::size_t, 3> pid3_type;
    typedef std::array<std::size_t, 4> pid4_type;
    typedef std::pair<pid2_type, potential_ptr> bond_potential_type;
    typedef std::pair<pid3_type, potential_ptr> angle_potential_type;
    typedef std::pair<pid4_type, potential_ptr> dihd_potential_type;
  
  public:

    LocalForceField() = default;
    ~LocalForceField() = default;

    void      calc_force(const particle_container_type& pcon);
    real_type calc_energy(const particle_container_type& pcon);

  private:

    std::vector<bond_potential_type>  bond_potentials;
    std::vector<angle_potential_type> angle_potentials;
    std::vector<dihd_potential_type>  dihd_potentials;
    BondLengthInteraction<traitsT>    bond;
    BondAngleInteraction<traitsT>     angle;
    DihedralAngleInteraction<traitsT> dihd;
};

template<typename traitsT>
void LocalForceField<traitsT>::calc_force(const particle_container_type& pcon)
{
    for(auto iter = bond_potentials.cbegin(); iter != bond_potentials.cend(); ++iter)
        bond.calc_force(pcon.at(iter->first[0]), pcon.at(iter->first[1]),
                        *(iter->second));

    for(auto iter = angle_potentials.cbegin(); iter != angle_potentials.cend(); ++iter)
        angle.calc_force(pcon.at(iter->first[0]), pcon.at(iter->first[1]),
                         pcon.at(iter->first[2]), *(iter->second));

    for(auto iter = dihd_potentials.cbegin(); iter != dihd_potentials.cend(); ++iter)
        dihd.calc_force(pcon.at(iter->first[0]), pcon.at(iter->first[1]),
                        pcon.at(iter->first[2]), pcon.at(iter->first[3]),
                        *(iter->second));

    return;
}

template<typename traitsT>
typename LocalForceField<traitsT>::real_type
LocalForceField<traitsT>::calc_energy(const particle_container_type& pcon)
{
    real_type energy = 0.0;
    for(auto iter = bond_potentials.cbegin(); iter != bond_potentials.cend(); ++iter)
        energy += bond.calc_energy(pcon.at(iter->first[0]), pcon.at(iter->first[1]),
                                   *(iter->second));

    for(auto iter = angle_potentials.cbegin(); iter != angle_potentials.cend(); ++iter)
        energy += angle.calc_energy(pcon.at(iter->first[0]), pcon.at(iter->first[1]),
                                    pcon.at(iter->first[2]), *(iter->second));

    for(auto iter = dihd_potentials.cbegin(); iter != dihd_potentials.cend(); ++iter)
        energy += dihd.calc_energy(pcon.at(iter->first[0]), pcon.at(iter->first[1]),
                                   pcon.at(iter->first[2]), pcon.at(iter->first[3]),
                                   *(iter->second));

    return energy;
}

} // mjolnir
#endif /* MJOLNIR_LOCAL_FORCE_FIELD */
