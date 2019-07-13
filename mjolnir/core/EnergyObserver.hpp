#ifndef MJOLNIR_CORE_ENERGY_OBSERVER_HPP
#define MJOLNIR_CORE_ENERGY_OBSERVER_HPP
#include <mjolnir/core/ObserverBase.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace mjolnir
{

template<typename traitsT>
class EnergyObserver final : public ObserverBase<traitsT>
{
  public:
    using base_type         = ObserverBase<traitsT>;
    using traits_type       = typename base_type::traits_type;
    using real_type         = typename base_type::real_type;
    using coordinate_type   = typename base_type::coordinate_type;
    using system_type       = typename base_type::system_type;
    using forcefield_type   = typename base_type::forcefield_type;

  public:

    explicit EnergyObserver(const std::string& filename_prefix)
      : prefix_(filename_prefix), file_name_(filename_prefix + ".ene")
    {
        // clear files and throw an error if the files cannot be opened.
        this->clear_file(this->file_name_);
    }
    ~EnergyObserver() override = default;

    void initialize(const std::size_t, const real_type,
                    const system_type&, const forcefield_type& ff) override
    {
        std::ofstream ofs(this->file_name_, std::ios::app);
        ofs << "# timestep  " << ff.list_energy_name() << " kinetic_energy\n";
        return;
    }

    void output(const std::size_t step, const real_type,
                const system_type& sys, const forcefield_type& ff) override
    {
        std::ofstream ofs(this->file_name_, std::ios::app);

        // if the width exceeds, operator<<(std::ostream, std::string) ignores
        // ostream::width and outputs whole string.
        ofs << std::setw(11) << std::left << std::to_string(step) << ' ';
        ofs << ff.dump_energy(sys) << ' ';
        ofs << std::setw(14) << std::right << this->calc_kinetic_energy(sys) << '\n';
        return;
    }

    void finalize(const std::size_t, const real_type,
                  const system_type&, const forcefield_type&) override
    {/* do nothing. */}

    std::string const& prefix() const noexcept override {return this->prefix_;}

  private:

    void clear_file(const std::string& fname) const
    {
        std::ofstream ofs(fname);
        if(not ofs.good())
        {
            throw_exception<std::runtime_error>(
                "[error] mjolnir::EnergyObserver: file open error: ", fname);
        }
        return;
    }

    real_type calc_kinetic_energy(const system_type& sys) const
    {
        real_type k = 0.0;
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            k += math::length_sq(sys.velocity(i)) * sys.mass(i);
        }
        return k * 0.5;
    }

  private:

    std::string prefix_;
    std::string file_name_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class EnergyObserver<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class EnergyObserver<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class EnergyObserver<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class EnergyObserver<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif


} // mjolnir
#endif//MJOLNIR_CORE_ENERGY_OBSERVER_HPP
