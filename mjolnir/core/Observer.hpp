#ifndef MJOLNIR_CORE_OBSERVER
#define MJOLNIR_CORE_OBSERVER
#include "Simulator.hpp"
#include <fstream>
#include <iomanip>

namespace mjolnir
{

template<typename traitsT>
class Observer
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;

  public:
    Observer() = default;
    Observer(const std::string& xyz, const std::string& ene)
        : xyz_name_(xyz), ene_name_(ene)
    {
        {
            std::ofstream ofs(xyz);
            if(not ofs.good())
                throw std::invalid_argument("file open error: " + xyz);
            ofs.close();
        }
        {
            std::ofstream ofs(ene);
            if(not ofs.good())
                throw std::invalid_argument("file open error: " + ene);
            ofs.close();
        }
    }
    ~Observer() = default;

    void output_coordinate(const Simulator<traitsT>& sim) const;
    void output_energy(const Simulator<traitsT>& sim) const;

    real_type calc_kinetic_energy(const ParticleContainer<traitsT>& pcon) const;

  private:

    std::string xyz_name_;
    std::string ene_name_;
};

template<typename traitsT>
inline void
Observer<traitsT>::output_coordinate(const Simulator<traitsT>& sim) const
{
    std::ofstream ofs(xyz_name_, std::ios::app);
    ofs << sim.particles().size() << std::endl;
    ofs << "time = " << sim.time() << std::endl;
    for(auto iter = sim.particles().cbegin(); iter != sim.particles().cend();
            ++iter)
    {
        ofs << "CA    " << std::fixed << std::setprecision(10)
            << iter->position[0] << " "
            << iter->position[1] << " "
            << iter->position[2] << " "
            << std::endl;
    }

    ofs.close();

    return ;
}

template<typename traitsT>
inline void
Observer<traitsT>::output_energy(const Simulator<traitsT>& sim) const
{
    std::ofstream ofs(ene_name_, std::ios::app);
    const real_type Ep = sim.calc_energy();
    const real_type Ek = this->calc_kinetic_energy(sim.particles());

    ofs << std::setw(8) << sim.time() << " "
        << std::setw(12) << Ep << " "
        << std::setw(12) << Ek << " "
        << std::setw(12) << Ep + Ek
        << std::endl;
    ofs.close();

    return ;
}

template<typename traitsT>
inline typename Observer<traitsT>::real_type
Observer<traitsT>::calc_kinetic_energy(const ParticleContainer<traitsT>& pcon) const
{
    real_type k = 0.0;
    for(auto iter = pcon.cbegin(); iter != pcon.cend(); ++iter)
        k += length_sq(iter->velocity) * iter->mass * 0.5;
    return k;
}



} // mjolnir
#endif /* MJOLNIR_CORE_OBSERVER */
