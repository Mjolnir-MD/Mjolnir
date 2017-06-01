#ifndef MJOLNIR_CORE_OBSERVER
#define MJOLNIR_CORE_OBSERVER
#include "System.hpp"
#include "ForceField.hpp"
#include <fstream>
#include <iomanip>

namespace mjolnir
{

template<typename traitsT>
class Observer
{
  public:
    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef ForceField<traits_type> forcefield_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:
    Observer(const std::string& xyz, const std::string& ene,
             const std::size_t interval)
        : interval_(interval), observe_count_(0),
          xyz_name_(xyz), ene_name_(ene)
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

    bool is_output_time()
    {
        ++observe_count_;
        if(this->observe_count_ == interval_)
        {
            this->observe_count_ = 0;
            return true;
        }
        return false;
    }

    void output(const real_type time, const system_type& sys,
                const forcefield_type& ff) const;

  private:

    real_type calc_kinetic_energy(const system_type& sys) const
    {
        real_type k = 0.0;
        for(const auto& particle : sys)
            k += length_sq(particle.velocity) * particle.mass;
        return k * 0.5;
    }

  private:

    std::size_t interval_;
    std::size_t observe_count_;
    std::string xyz_name_;
    std::string ene_name_;
};

template<typename traitsT>
inline void Observer<traitsT>::output(
    const real_type time, const system_type& sys, const forcefield_type& ff) const
{
    std::ofstream ofs(xyz_name_, std::ios::app);
    ofs << sys.size() << "\n" << time << "\n";
    for(const auto& particle : sys)
    {// TODO change output format
        ofs << "CA    " << std::fixed << std::setprecision(8)
            << particle.position[0] << " "
            << particle.position[1] << " "
            << particle.position[2] << std::endl;
    }
    ofs.close();

    // TODO separate energy terms
    ofs.open(ene_name_, std::ios::app);
    ofs << time << " " << ff.calc_energy(sys) << std::endl;
    ofs.close();

    return ;
}

} // mjolnir
#endif /* MJOLNIR_CORE_OBSERVER */
