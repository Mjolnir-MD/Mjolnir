#ifndef MJOLNIR_CORE_TRR_OBSERVER_HPP
#define MJOLNIR_CORE_TRR_OBSERVER_HPP
#include <mjolnir/core/ObserverBase.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/Unit.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

// It is an observer class that observes the current state of system and outputs
// the status to a file.
// TRRObserver outputs .trr file.
//
// XXX Here the unit system is defined as the current simulator setup.
//     No unit conversions are applied.

namespace mjolnir
{

template<typename traitsT>
class TRRObserver final : public ObserverBase<traitsT>
{
  public:
    using base_type         = ObserverBase<traitsT>;
    using traits_type       = typename base_type::traits_type;
    using real_type         = typename base_type::real_type;
    using coordinate_type   = typename base_type::coordinate_type;
    using system_type       = typename base_type::system_type;
    using forcefield_type   = typename base_type::forcefield_type;

  public:

    explicit TRRObserver(const std::string& filename_prefix)
      : base_type(), prefix_(filename_prefix),
        trr_name_(filename_prefix + std::string(".trr"))
    {
        // clear files and throw an error if the files cannot be opened.
        this->clear_file(this->trr_name_);
    }
    ~TRRObserver() override = default;

    void initialize(const std::size_t total_step, const real_type dt,
                    const system_type& sys, const forcefield_type& ff) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        if(physics::constants<real_type>::m_to_length() !=
              unit::constants<real_type>::m_to_nm)
        {
            MJOLNIR_LOG_NOTICE("length unit seems to be different from [nm]."    );
            MJOLNIR_LOG_NOTICE("The resulting .trr file might look weird because");
            MJOLNIR_LOG_NOTICE("of the unit conversion (normally [nm] is used)"  );
        }
        return; // do nothing.
    }

    void output(const std::size_t step, const real_type dt,
                const system_type& sys, const forcefield_type& ff) override;

    void finalize(const std::size_t, const real_type dt,
                  const system_type&, const forcefield_type&) override
    {
        return; // do nothing.
    }

    std::string const& prefix() const noexcept override {return prefix_;}

  private:

    void clear_file(const std::string& fname) const
    {
        std::ofstream ofs(fname);
        if(not ofs.good())
        {
            throw_exception<std::runtime_error>("[error] mjolnir::TRRObserver: "
                    "file open error: ", fname);
        }
        return;
    }

    // it is a helper function to write unitcell region size.
    // for UnlimitedBoundary, returns zero.
    static std::int32_t boxinfo_size(
        const UnlimitedBoundary<real_type, coordinate_type>&) noexcept
    {
        // no unitcell information needed. block size is none.
        return 0;
    }
    // for CuboidalPeriodicBoundary, returns .
    static std::int32_t boxinfo_size(
        const CuboidalPeriodicBoundary<real_type, coordinate_type>&) noexcept
    {
        // unitcell information required.
        // write 3 * 3D vectors for each edge vector of the cuboid.
        return sizeof(real_type) * 3 * 3;
    }

    // it is a helper function to write unitcell region size.
    // for UnlimitedBoundary, returns zero.
    static void write_boxinfo(std::ostream&,
        const UnlimitedBoundary<real_type, coordinate_type>&) noexcept
    {
        // no unitcell information needed. block size is none.
        return ;
    }
    // for CuboidalPeriodicBoundary, returns .
    static void write_boxinfo(std::ostream& os,
        const CuboidalPeriodicBoundary<real_type, coordinate_type>& bdry) noexcept
    {
        const coordinate_type& width = bdry.width();

        // vector for x direction
        detail::write_as_bytes(os, real_type(math::X(width)));
        detail::write_as_bytes(os, real_type(0));
        detail::write_as_bytes(os, real_type(0));

        // vector for y direction
        detail::write_as_bytes(os, real_type(0));
        detail::write_as_bytes(os, real_type(math::Y(width)));
        detail::write_as_bytes(os, real_type(0));

        // vector for z direction
        detail::write_as_bytes(os, real_type(0));
        detail::write_as_bytes(os, real_type(0));
        detail::write_as_bytes(os, real_type(math::Z(width)));
        return ;
    }

  private:

    std::string prefix_;
    std::string trr_name_;
    std::size_t number_of_frames_;
};

template<typename traitsT>
inline void TRRObserver<traitsT>::output(
    const std::size_t step, const real_type dt,
    const system_type& sys, const forcefield_type& ff)
{
    using self_type = TRRObserver<traitsT>;
    std::ofstream ofs(this->trr_name_, std::ios::app | std::ios::binary);

    // ------------------------------------------------------------------------
    // write frame header
    detail::write_as_bytes(ofs, std::int32_t(1993)); // magic number for trr
    detail::write_as_bytes(ofs, std::int32_t(  13));

    detail::write_as_bytes(ofs, std::int32_t(this->prefix_.size()));
    ofs.write(this->prefix_.data(), this->prefix_.size()); // title

    const std::int32_t sys_size = sys.size(); // number of particles
    const std::int32_t crd_size = sizeof(real_type) * 3 * sys_size;

    detail::write_as_bytes(ofs, std::int32_t(0));
    detail::write_as_bytes(ofs, std::int32_t(0));
    detail::write_as_bytes(ofs, self_type::boxinfo_size(sys.boundary()));
    detail::write_as_bytes(ofs, std::int32_t(0));
    detail::write_as_bytes(ofs, std::int32_t(0));
    detail::write_as_bytes(ofs, std::int32_t(0));
    detail::write_as_bytes(ofs, std::int32_t(0));
    detail::write_as_bytes(ofs, crd_size       );    // xyz size
    detail::write_as_bytes(ofs, crd_size       );    // vel size
    detail::write_as_bytes(ofs, crd_size       );    // force size
    detail::write_as_bytes(ofs, sys_size       );    // natoms
    detail::write_as_bytes(ofs, std::int32_t(step)); // step
    detail::write_as_bytes(ofs, std::int32_t(0));

    detail::write_as_bytes(ofs,             dt); // tstep
    detail::write_as_bytes(ofs, real_type(0.0));

    // ------------------------------------------------------------------------
    // write box size

    self_type::write_boxinfo(ofs, sys.boundary());

    // ------------------------------------------------------------------------
    // write positions

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        detail::write_as_bytes(ofs, math::X(sys.position(i)));
        detail::write_as_bytes(ofs, math::Y(sys.position(i)));
        detail::write_as_bytes(ofs, math::Z(sys.position(i)));
    }

    // ------------------------------------------------------------------------
    // write velocities
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        detail::write_as_bytes(ofs, math::X(sys.velocity(i)));
        detail::write_as_bytes(ofs, math::Y(sys.velocity(i)));
        detail::write_as_bytes(ofs, math::Z(sys.velocity(i)));
    }

    // ------------------------------------------------------------------------
    // write forces
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        detail::write_as_bytes(ofs, math::X(sys.force(i)));
        detail::write_as_bytes(ofs, math::Y(sys.force(i)));
        detail::write_as_bytes(ofs, math::Z(sys.force(i)));
    }

    return ;
}

} // mjolnir
#endif // MJOLNIR_CORE_TRR_OBSERVER_HPP
