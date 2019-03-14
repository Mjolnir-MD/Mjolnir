#ifndef MJOLNIR_CORE_XYZ_OBSERVER_HPP
#define MJOLNIR_CORE_XYZ_OBSERVER_HPP
#include <mjolnir/core/ObserverBase.hpp>
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

  public:

    explicit XYZObserver(const std::string& filename_prefix)
      : base_type(), prefix_(filename_prefix),
        xyz_name_(filename_prefix + std::string("_position.xyz")),
        vel_name_(filename_prefix + std::string("_velocity.xyz")),
    {
        // clear files and throw an error if the files cannot be opened.
        this->clear_file(this->xyz_name_);
        this->clear_file(this->vel_name_);
    }
    ~XYZObserver() override = default;

    void initialize(const std::size_t total_step,
                    const system_type& sys, const forcefield_type& ff) override
    {
        // do nothing.
    }

    void output(const std::size_t step,
                const system_type& sys, const forcefield_type& ff) override
    {
        // -------------------------------------------------------------------
        // output positions
        {
            std::ofstream ofs(xyz_name_, std::ios::app);
            ofs << sys.size() << "\nstep = " << step << '\n';
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                const auto& p = sys.position(i);
                ofs << sys.name(i) << ' ' << std::fixed << std::setprecision(8)
                    << math::X(p)  << ' ' << math::Y(p) << ' ' << math::Z(p)
                    << '\n';
            }
            ofs.close();
        }

        // -------------------------------------------------------------------
        // output velocities
        {
            std::ofstream ofs(vel_name_, std::ios::app);
            ofs << sys.size() << "\nstep = " << step << '\n';
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                const auto& v = sys.velocity(i);
                ofs << sys.name(i) << ' ' << std::fixed << std::setprecision(8)
                    << math::X(v)  << ' ' << math::Y(v) << ' ' << math::Z(v)
                    << '\n';
            }
            ofs.close();
        }
        return ;
    }

    void finalize(const std::size_t total_step,
                  const system_type& sys, const forcefield_type& ff) override
    {
        // do nothing.
    }

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

  private:

    std::string prefix_;
    std::string xyz_name_;
    std::string vel_name_;
};

} // mjolnir
#endif // MJOLNIR_CORE_XYZ_OBSERVER_HPP
