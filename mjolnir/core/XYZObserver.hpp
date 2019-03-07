#ifndef MJOLNIR_CORE_XYZ_OBSERVER_HPP
#define MJOLNIR_CORE_XYZ_OBSERVER_HPP
#include <mjolnir/core/ObserverBase.hpp>
#include <mjolnir/util/progress_bar.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace mjolnir
{

template<typename traitsT>
class XYZObserver final : public ObserverBase<traitsT>
{
  public:
    using base_type         = ObserverBase<traitsT>;
    using traits_type       = typename base_type::traits_type;
    using real_type         = typename base_type::real_type;
    using coordinate_type   = typename base_type::coordinate_type;
    using system_type       = typename base_type::system_type;
    using forcefield_type   = typename base_type::forcefield_type;
    using progress_bar_type = progress_bar<50>;

  public:

    XYZObserver(const std::string& filename_prefix, bool output_progress = false)
      : base_type(), output_progress_(output_progress), progress_bar_(1),
        prefix_(filename_prefix),
        xyz_name_(filename_prefix + std::string("_position.xyz")),
        vel_name_(filename_prefix + std::string("_velocity.xyz")),
        ene_name_(filename_prefix + std::string(".ene"))
    {
        // clear files and throw an error if the files cannot be opened.
        this->clear_file(this->xyz_name_);
        this->clear_file(this->vel_name_);
        this->clear_file(this->ene_name_);
    }
    ~XYZObserver() override = default;

    void initialize(const std::size_t total_step,
                    const system_type& sys, const forcefield_type& ff) override
    {
        this->progress_bar_.reset(total_step); // set total_step

        std::ofstream ofs(this->ene_name_, std::ios::app);
        ofs << "# timestep  " << ff.list_energy_name() << " kinetic_energy\n";
        return;
    }

    void output(const std::size_t step,
                const system_type& sys, const forcefield_type& ff) override;

    std::string const& prefix() const noexcept override {return prefix_;}

  private:

    void clear_file(const std::string& fname) const
    {
        std::ofstream ofs(fname);
        if(not ofs.good())
        {
            throw_exception<std::runtime_error>("[error] mjolnir::XYZObserver: "
                    "file open error: ", fname);
        }
        return;
    }

    real_type calc_kinetic_energy(const system_type& sys) const
    {
        real_type k = 0.0;
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            k += math::length_sq(sys[i].velocity) * sys[i].mass;
        }
        return k * 0.5;
    }

  private:

    bool output_progress_;
    std::string prefix_;
    std::string xyz_name_;
    std::string vel_name_;
    std::string ene_name_;
    progress_bar_type progress_bar_;
};

template<typename traitsT>
inline void XYZObserver<traitsT>::output(
    const std::size_t step, const system_type& sys, const forcefield_type& ff)
{
    std::ofstream ofs(xyz_name_, std::ios::app);
    ofs << sys.size() << "\nstep = " << step << '\n';
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto& p = sys.position(i);
        ofs << sys.name(i) << ' ' << std::fixed << std::setprecision(8)
            << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';
    }
    ofs.close();

    ofs.open(vel_name_, std::ios::app);
    ofs << sys.size() << "\nstep = " << step << '\n';
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto& v = sys.velocity(i);
        ofs << sys.name(i) << ' ' << std::fixed << std::setprecision(8)
            << v[0] << ' ' << v[1] << ' ' << v[2] << '\n';
    }
    ofs.close();

    ofs.open(ene_name_, std::ios::app);
    // if the width exceeds, operator<<(std::ostream, std::string) ignores
    // ostream::width and outputs whole string.
    ofs << std::setw(11) << std::left << std::to_string(step) << ' ';
    ofs << ff.dump_energy(sys) << ' ';
    ofs << std::setw(14) << std::right << this->calc_kinetic_energy(sys) << '\n';
    ofs.close();

    // XXX consider introducing template argument to remove this if-branching
    //     at the compile time
    if(this->output_progress_)
    {
        std::cerr << progress_bar_.format(step);
        if(step == progress_bar_.total()){std::cerr << std::endl;}
    }
    return ;
}

} // mjolnir
#endif // MJOLNIR_CORE_XYZ_OBSERVER_HPP
