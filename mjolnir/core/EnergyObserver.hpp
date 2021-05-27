#ifndef MJOLNIR_CORE_ENERGY_OBSERVER_HPP
#define MJOLNIR_CORE_ENERGY_OBSERVER_HPP
#include <mjolnir/util/is_finite.hpp>
#include <mjolnir/core/ObserverBase.hpp>
#include <mjolnir/core/Unit.hpp>
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
    using boundary_type     = typename traits_type::boundary_type;

  public:

    explicit EnergyObserver(const std::string& filename_prefix)
      : prefix_(filename_prefix), file_name_(filename_prefix + ".ene")
    {
        // clear files and throw an error if the files cannot be opened.
        this->clear_file(this->file_name_);
    }
    ~EnergyObserver() override {}

    void initialize(const std::size_t, const std::size_t, const real_type,
                    const system_type& sys, const forcefield_type& ff) override
    {
        using phys_constants = physics::constants<real_type>;
        std::ofstream ofs(this->file_name_, std::ios::app);
        ofs << "# unit of length : " << phys_constants::length_unit()
            << ", unit of energy : " << phys_constants::energy_unit() << '\n';
        ofs << "# timestep  ";

        std::string names;
        ff->format_energy_name(names); // write into it
        ofs << names;

        ofs << " kinetic_energy             Pxx             Pyy             Pzz";
        for(const auto& attr : sys.attributes())
        {
            ofs << " attribute:" << attr.first;
        }
        ofs << '\n';
        return;
    }

    // update column names and widths if forcefield changed.
    void update(const std::size_t,  const real_type,
                const system_type& sys, const forcefield_type& ff) override
    {
        std::ofstream ofs(this->file_name_, std::ios::app);
        ofs << "# timestep  ";

        std::string names;
        ff->format_energy_name(names); // write into it
        ofs << names;

        ofs << " kinetic_energy             Pxx             Pyy             Pzz";
        for(const auto& attr : sys.attributes())
        {
            ofs << " attribute:" << attr.first;
        }
        ofs << '\n';
        return;
    }

    void output(const std::size_t step, const real_type,
                const system_type& sys, const forcefield_type& ff) override
    {
        bool is_ok = true;
        std::ofstream ofs(this->file_name_, std::ios::app);

        // if the width exceeds, operator<<(std::ostream, std::string) ignores
        // ostream::width and outputs whole string.
        ofs << std::setw(11) << std::left << std::to_string(step) << ' ';

        std::string energies;
        const auto potential_energy = ff->format_energy(sys, energies);
        ofs << energies;
        if(!is_finite(potential_energy))
        {
            MJOLNIR_GET_DEFAULT_LOGGER();
            MJOLNIR_LOG_ERROR("potential energy becomes NaN.");
            is_ok = false;
        }

        real_type Ek(0), Px(0), Py(0), Pz(0);
        std::tie(Ek, Px, Py, Pz) = this->calc_energy_and_pressure(sys);

        ofs << std::setw(14) << std::right << std::fixed << Ek;
        ofs << std::setw(16) << std::right << std::setprecision(8) << std::scientific << Px;
        ofs << std::setw(16) << std::right << std::setprecision(8) << std::scientific << Py;
        ofs << std::setw(16) << std::right << std::setprecision(8) << std::scientific << Pz;
        if(!is_finite(Ek))
        {
            MJOLNIR_GET_DEFAULT_LOGGER();
            MJOLNIR_LOG_ERROR("kinetic energy becomes NaN.");
            is_ok = false;
        }

        for(const auto& attr : sys.attributes())
        {
            if(!is_finite(attr.second))
            {
                MJOLNIR_GET_DEFAULT_LOGGER();
                MJOLNIR_LOG_ERROR(attr.first, " becomes NaN.");
                is_ok = false;
            }
            ofs << ' ' << std::setw(10 + attr.first.size()) << std::right
                << std::fixed << attr.second;
        }
        ofs << std::endl; // flush before throwing an exception

        if(!is_ok)
        {
            throw std::runtime_error("Energy value becomes NaN");
        }
        return;
    }

    void finalize(const std::size_t, const real_type,
                  const system_type&, const forcefield_type&) override
    {/* do nothing. */}

    std::string const& prefix() const noexcept override {return this->prefix_;}

  private:

    std::tuple<real_type, real_type, real_type, real_type>
    calc_energy_and_pressure(const system_type& sys)
    {
        const auto cell_width = sys.boundary().width();
        const auto volume = math::X(cell_width) * math::Y(cell_width) * math::Z(cell_width);
        const auto rvolume = real_type(1) / volume;

        real_type Ek(0);
        real_type Px(sys.virial()(0, 0));
        real_type Py(sys.virial()(1, 1));
        real_type Pz(sys.virial()(2, 2));

        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto  m = sys.mass(i);
            const auto& v = sys.velocity(i);

            Ek += math::length_sq(v) * m;
            Px += m * math::X(v) * math::X(v);
            Py += m * math::Y(v) * math::Y(v);
            Pz += m * math::Z(v) * math::Z(v);
        }
        Ek *= 0.5;
        Px *= rvolume;
        Py *= rvolume;
        Pz *= rvolume;
        return std::make_tuple(Ek, Px, Py, Pz);
    }

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

  private:

    std::string prefix_;
    std::string file_name_;
    std::vector<std::size_t> widths_; // column width to format energy values
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class EnergyObserver<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class EnergyObserver<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class EnergyObserver<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class EnergyObserver<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif


} // mjolnir
#endif//MJOLNIR_CORE_ENERGY_OBSERVER_HPP
