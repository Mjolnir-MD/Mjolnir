#ifndef MJOLNIR_CORE_OBSERVER
#define MJOLNIR_CORE_OBSERVER
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceField.hpp>
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
    Observer(const std::string& filename_prefix) :
        xyz_name_(filename_prefix + std::string(".xyz")),
        ene_name_(filename_prefix + std::string(".ene"))
    {
        // clear the contents
        {
            std::ofstream ofs(this->xyz_name_);
            if(not ofs.good())
            {
                throw std::runtime_error("file open error: " + this->xyz_name_);
            }
            ofs.close();
        }
        {
            std::ofstream ofs(this->ene_name_);
            if(not ofs.good())
            {
                throw std::runtime_error("file open error: " + this->ene_name_);
            }
            ofs.close();
        }
    }
    ~Observer() = default;

    void initialize(const system_type& sys, const forcefield_type& ff) const
    {
        std::ofstream ofs(this->ene_name_, std::ios::app);
        ofs << '#' << ff.list_energy_name() << "kinetic_energy\n";
        return;
    }

    void output(const real_type time, const system_type& sys,
                const forcefield_type& ff) const;

  private:

    real_type calc_kinetic_energy(const system_type& sys) const
    {
        real_type k = 0.0;
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            k += length_sq(sys[i].velocity) * sys[i].mass;
        }
        return k * 0.5;
    }

  private:

    std::string xyz_name_;
    std::string ene_name_;
};

template<typename traitsT>
inline void Observer<traitsT>::output(
    const real_type time, const system_type& sys, const forcefield_type& ff) const
{
    std::ofstream ofs(xyz_name_, std::ios::app);
    ofs << sys.size() << '\n' << time << '\n';
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto& p = sys.position(i);
        ofs << "CA    " << std::fixed << std::setprecision(8)
            << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
    }
    ofs.close();

    ofs.open(ene_name_, std::ios::app);
    ofs << time << ' ' << ff.dump_energy(sys) << ' '
        << this->calc_kinetic_energy(sys) << '\n';
    ofs.close();
    return ;
}

} // mjolnir
#endif /* MJOLNIR_CORE_OBSERVER */
