#ifndef MJOLNIR_CORE_DCD_OBSERVER_HPP
#define MJOLNIR_CORE_DCD_OBSERVER_HPP
#include <mjolnir/core/ObserverBase.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/System.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>

namespace mjolnir
{

template<typename traitsT>
class DCDObserver final : public ObserverBase<traitsT>
{
  public:
    using base_type         = ObserverBase<traitsT>;
    using traits_type       = typename base_type::traits_type;
    using real_type         = typename base_type::real_type;
    using coordinate_type   = typename base_type::coordinate_type;
    using system_type       = typename base_type::system_type;
    using forcefield_type   = typename base_type::forcefield_type;

  public:

    explicit DCDObserver(const std::string& filename_prefix)
      : base_type(), prefix_(filename_prefix),
        pos_name_(filename_prefix + std::string("_position.dcd")),
        vel_name_(filename_prefix + std::string("_velocity.dcd"))
    {
        // clear files and throw an error if the files cannot be opened.
        this->clear_file(this->pos_name_);
        this->clear_file(this->vel_name_);
    }
    ~DCDObserver() override = default;

    void initialize(const std::size_t total_step, const real_type dt,
                    const system_type& sys, const forcefield_type& ff) override
    {
        this->write_header(this->pos_name_, total_step, dt, sys, ff);
        this->write_header(this->vel_name_, total_step, dt, sys, ff);

        // buffer to convert sys and dcd format
        this->buffer_x_.resize(sys.size());
        this->buffer_y_.resize(sys.size());
        this->buffer_z_.resize(sys.size());
        return;
    }

    void update(const std::size_t,      const real_type,
                const system_type& sys, const forcefield_type&) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        if(sys.size() != this->buffer_x_.size())
        {
            MJOLNIR_LOG_NOTICE("the number of particles changed.");
            MJOLNIR_LOG_NOTICE("Most of dcd file readers assumes it is a constant.");
            MJOLNIR_LOG_NOTICE("It may cause some problems.");

            this->buffer_x_.resize(sys.size());
            this->buffer_y_.resize(sys.size());
            this->buffer_z_.resize(sys.size());
        }
        return;
    }

    void output(const std::size_t step, const real_type dt,
                const system_type& sys, const forcefield_type& ff) override;

    void finalize(const std::size_t, const real_type,
                  const system_type&, const forcefield_type&) override
    {
        // update # of frames in the header region
        {
            // both `in` and `out` flags are needed to keep the other parts of
            // file. otherwise, we will lost almost everything by overwriting.
            // Also `ate` should be used instead of `app`, because `app` moves
            // the position automatically at the end when we call `write`.

            std::ofstream ofs(this->pos_name_,
                std::ios::binary | std::ios::ate | std::ios::in | std::ios::out);
            // skip the first block size and the signature "CORD"
            ofs.seekp(2 * sizeof(std::int32_t), std::ios::beg);

            const std::int32_t number_of_frames(this->number_of_frames_);
            detail::write_as_bytes(ofs, number_of_frames);
        }
        {
            std::ofstream ofs(this->vel_name_,
                std::ios::binary | std::ios::ate | std::ios::in | std::ios::out);
            // skip the first block size and the signature "CORD"
            ofs.seekp(2 * sizeof(std::int32_t), std::ios::beg);

            const std::int32_t number_of_frames(this->number_of_frames_);
            detail::write_as_bytes(ofs, number_of_frames);
        }
        return;
    }

    std::string const& prefix() const noexcept override {return prefix_;}

  private:

    void clear_file(const std::string& fname) const
    {
        std::ofstream ofs(fname);
        if(not ofs.good())
        {
            throw_exception<std::runtime_error>("[error] mjolnir::DCDObserver: "
                    "file open error: ", fname);
        }
        return;
    }

    void write_header(const std::string& fname,
                      const std::size_t  total_step_sz, const real_type dt,
                      const system_type& sys, const forcefield_type&) const
    {
        std::ofstream ofs(fname, std::ios::binary | std::ios::app);
        if(not ofs.good())
        {
            throw_exception<std::runtime_error>(
                    "[error] mjolnir::DCDObserver: file open error: ", fname);
        }

        /* the first block */
        {
            const std::int32_t block_size(84);
            detail::write_as_bytes(ofs, block_size);
            ofs.write("CORD", 4);

            const std::int32_t number_of_frames(0);
            detail::write_as_bytes(ofs, number_of_frames);

            const std::int32_t index_of_first(0);
            detail::write_as_bytes(ofs, index_of_first);

            const std::int32_t save_interval(0);
            detail::write_as_bytes(ofs, save_interval);

            const std::int32_t total_step(total_step_sz);
            detail::write_as_bytes(ofs, total_step);

            const std::int32_t total_chains(sys.topology().number_of_molecules());
            detail::write_as_bytes(ofs, total_chains);

            const std::int32_t zero(0);
            // 4 * integers with null flag
            for(std::size_t i=0; i<4; ++i) {detail::write_as_bytes(ofs, zero);}

            const float delta_t(dt);
            detail::write_as_bytes(ofs, delta_t);

            const std::int32_t has_unitcell =
                DCDObserver<traitsT>::unitcell_flag(sys.boundary());
            detail::write_as_bytes(ofs, has_unitcell);

            // 8 * integers with null flag
            for(std::size_t i=0; i<8; ++i) {detail::write_as_bytes(ofs, zero);}

            const std::int32_t version(24);
            detail::write_as_bytes(ofs, version);

            detail::write_as_bytes(ofs, block_size);
        }

        /* the second block */
        {
            const std::int32_t block_size(84);
            detail::write_as_bytes(ofs, block_size);

            const std::int32_t number_of_lines(1);
            detail::write_as_bytes(ofs, number_of_lines);

            const char comment[80] = "Mjolnir -- copyright (c) Toru Niina 2016"
                                     "-now distributed under the MIT License.";
            ofs.write(comment, 80);

            detail::write_as_bytes(ofs, block_size);
        }

        /* the third block */
        {
            const std::int32_t block_size(4);
            detail::write_as_bytes(ofs, block_size);

            const std::int32_t number_of_particles(sys.size());
            detail::write_as_bytes(ofs, number_of_particles);

            detail::write_as_bytes(ofs, block_size);
        }
        return;
    }

    // it is a helper function to write unitcell flags.
    // for UnlimitedBoundary, returns zero.
    static std::int32_t unitcell_flag(
        const UnlimitedBoundary<real_type, coordinate_type>&) noexcept
    {
        // no unitcell information needed. disable the flag
        return 0;
    }
    // for CuboidalPeriodicBoundary, returns one.
    static std::int32_t unitcell_flag(
        const CuboidalPeriodicBoundary<real_type, coordinate_type>&) noexcept
    {
        // unitcell information required. turn the flag on.
        return 1;
    }

    // it is a helper function to write unitcell block if needed.
    // for UnlimitedBoundary, do nothing.
    static void write_unitcell_if_needed(std::ostream&,
        const UnlimitedBoundary<real_type, coordinate_type>&) noexcept
    {
        return ; // do nothing. boundary does not exists.
    }
    // for CuboidalPeriodicBoundary, writes the boundary width and angles
    static void write_unitcell_if_needed(std::ostream& os,
        const CuboidalPeriodicBoundary<real_type, coordinate_type>& boundary) noexcept
    {
        // unit cell length
        const double A = math::X(boundary.width());
        const double B = math::Y(boundary.width());
        const double C = math::Z(boundary.width());

        // angles are always 90 degree because it's cuboid.
        // for earlier versions, it was a cosine value of the angle.
        // Now it accepts degrees. For the clarity, I use degrees here.
        const double alpha = 90.0;
        const double beta  = 90.0;
        const double gamma = 90.0;

        const std::int32_t block_size = sizeof(double) * 6;
        detail::write_as_bytes(os, block_size);

        // I'm serious. the order is correct.
        detail::write_as_bytes(os, A    );
        detail::write_as_bytes(os, gamma);
        detail::write_as_bytes(os, B    );
        detail::write_as_bytes(os, beta );
        detail::write_as_bytes(os, alpha);
        detail::write_as_bytes(os, C    );

        detail::write_as_bytes(os, block_size);
        return ;
    }

  private:

    std::string prefix_;
    std::string pos_name_;
    std::string vel_name_;
    std::size_t number_of_frames_;
    std::vector<float> buffer_x_;
    std::vector<float> buffer_y_;
    std::vector<float> buffer_z_;
};

template<typename traitsT>
inline void DCDObserver<traitsT>::output(
    const std::size_t, const real_type,
    const system_type& sys, const forcefield_type&)
{
    number_of_frames_ += 1;
    assert(this->buffer_x_.size() == sys.size());
    assert(this->buffer_y_.size() == sys.size());
    assert(this->buffer_z_.size() == sys.size());

    // ------------------------------------------------------------------------
    // write position
    {
        std::ofstream ofs(this->pos_name_, std::ios::app | std::ios::binary);

        DCDObserver<traitsT>::write_unitcell_if_needed(ofs, sys.boundary());

        for(std::size_t i=0; i<sys.size(); ++i)
        {
            this->buffer_x_[i] = static_cast<float>(math::X(sys.position(i)));
            this->buffer_y_[i] = static_cast<float>(math::Y(sys.position(i)));
            this->buffer_z_[i] = static_cast<float>(math::Z(sys.position(i)));
        }
        const std::int32_t block_size(sizeof(float) * sys.size());
        {
            detail::write_as_bytes(ofs, block_size);
            ofs.write(reinterpret_cast<const char*>(this->buffer_x_.data()),
                      block_size);
            detail::write_as_bytes(ofs, block_size);
        }
        {
            detail::write_as_bytes(ofs, block_size);
            ofs.write(reinterpret_cast<const char*>(this->buffer_y_.data()),
                      block_size);
            detail::write_as_bytes(ofs, block_size);
        }
        {
            detail::write_as_bytes(ofs, block_size);
            ofs.write(reinterpret_cast<const char*>(this->buffer_z_.data()),
                      block_size);
            detail::write_as_bytes(ofs, block_size);
        }
    }

    // ------------------------------------------------------------------------
    // write velocity
    {
        std::ofstream ofs(this->vel_name_, std::ios::app | std::ios::binary);
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            this->buffer_x_[i] = static_cast<float>(math::X(sys.velocity(i)));
            this->buffer_y_[i] = static_cast<float>(math::Y(sys.velocity(i)));
            this->buffer_z_[i] = static_cast<float>(math::Z(sys.velocity(i)));
        }
        const std::int32_t block_size(sizeof(float) * sys.size());
        {
            detail::write_as_bytes(ofs, block_size);
            ofs.write(reinterpret_cast<const char*>(this->buffer_x_.data()),
                      block_size);
            detail::write_as_bytes(ofs, block_size);
        }
        {
            detail::write_as_bytes(ofs, block_size);
            ofs.write(reinterpret_cast<const char*>(this->buffer_y_.data()),
                      block_size);
            detail::write_as_bytes(ofs, block_size);
        }
        {
            detail::write_as_bytes(ofs, block_size);
            ofs.write(reinterpret_cast<const char*>(this->buffer_z_.data()),
                      block_size);
            detail::write_as_bytes(ofs, block_size);
        }
    }
    return ;
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class DCDObserver<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class DCDObserver<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class DCDObserver<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class DCDObserver<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif // MJOLNIR_CORE_DCD_OBSERVER_HPP
