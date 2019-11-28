#ifndef MJOLNIR_CORE_TRR_LOADER_HPP
#define MJOLNIR_CORE_TRR_LOADER_HPP
#include <mjolnir/core/LoaderBase.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/util/logger.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace mjolnir
{

template<typename traitsT>
class TRRLoader final : public LoaderBase<traitsT>
{
  public:
    using base_type         = LoaderBase<traitsT>;
    using traits_type       = typename base_type::traits_type;
    using real_type         = typename base_type::real_type;
    using coordinate_type   = typename base_type::coordinate_type;
    using system_type       = typename base_type::system_type;

  public:

    explicit TRRLoader(const std::string& filename)
        : base_type(), has_unitcell_(false), filename_(filename),
          number_of_frames_(0), number_of_particles_(0),
          position_block_size_(0), velocity_block_size_(0), force_block_size_(0),
          file_(filename_, std::ios::binary | std::ios::in)
    {
        if(!file_.good())
        {
            throw_exception<std::runtime_error>("[error] mjolnir::TRRLoader: "
                    "file open error: ", filename_);
        }
    }
    ~TRRLoader() override {}

    void initialize() override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        if(physics::constants<real_type>::m_to_length() !=
              unit::constants<real_type>::m_to_nm())
        {
            MJOLNIR_LOG_WARN("length unit seems to be different from [nm].");
            MJOLNIR_LOG_WARN("But it does not do any unit conversion while "
                             "reading it because TRRObserver does so.");
        }

        // -------------------------------------------------------------------
        // calculate file size to estimate total number of steps

        file_.seekg(0, std::ios_base::end);
        const auto at_end = file_.tellg();
        file_.seekg(0, std::ios_base::beg);
        const auto file_size = at_end - file_.tellg();

        this->file_.clear(); // clear bitflags
        this->file_.peek();  // update bitflags

        // -------------------------------------------------------------------
        // read the first frame header and store information.

        std::ignore = detail::read_bytes_as<std::int32_t>(file_); // magic number
        std::ignore = detail::read_bytes_as<std::int32_t>(file_); // version
        const auto len_title = detail::read_bytes_as<std::int32_t>(file_);
        std::vector<char> title(len_title);
        file_.read(title.data(), len_title);
        {
            std::string t(title.begin(), title.end());
            MJOLNIR_LOG_NOTICE("TRR file title is ", t);
        }
        const auto   ir_size = detail::read_bytes_as<std::int32_t>(file_);
        const auto    e_size = detail::read_bytes_as<std::int32_t>(file_);
        const auto  box_size = detail::read_bytes_as<std::int32_t>(file_);
        const auto  vir_size = detail::read_bytes_as<std::int32_t>(file_);
        const auto pres_size = detail::read_bytes_as<std::int32_t>(file_);
        const auto  top_size = detail::read_bytes_as<std::int32_t>(file_);
        const auto  sym_size = detail::read_bytes_as<std::int32_t>(file_);
        this->position_block_size_ = detail::read_bytes_as<std::int32_t>(file_);
        this->velocity_block_size_ = detail::read_bytes_as<std::int32_t>(file_);
        this->force_block_size_    = detail::read_bytes_as<std::int32_t>(file_);
        this->number_of_particles_ = detail::read_bytes_as<std::int32_t>(file_);
        std::ignore = detail::read_bytes_as<std::int32_t>(file_); // step
        std::ignore = detail::read_bytes_as<std::int32_t>(file_); // nre
        std::ignore = detail::read_bytes_as<real_type>(file_);    // t
        std::ignore = detail::read_bytes_as<real_type>(file_);    // lambda

        // TODO: check floating point precision before reading real type
        //       currently it assumes that the precision is correctly defined

        // -------------------------------------------------------------------
        // estimate total trajectory length

        const auto snapshot_size = 64 + len_title + 2 * sizeof(real_type) +
            ir_size + e_size + box_size + vir_size + pres_size + top_size +
            sym_size + position_block_size_ + velocity_block_size_ +
            force_block_size_;

        this->number_of_frames_ = file_size / snapshot_size;

        if(file_size % snapshot_size != 0)
        {
            MJOLNIR_LOG_WARN("File length is not a multiple of the size of the"
                             " first snapshot!");
        }
        MJOLNIR_LOG_NOTICE("Estimated number of frames in ", filename_, " is ",
                           this->number_of_frames_);

        // -------------------------------------------------------------------
        // rewind and reset bitflags inside

        this->file_.clear();                      // clear bitflags
        this->file_.seekg(0, std::ios_base::beg); // rewind
        this->file_.peek();                       // update bitflags
        return;
    }

    std::size_t num_particles() const noexcept {return number_of_particles_;}
    std::size_t num_frames()    const noexcept {return number_of_frames_;}
    bool        is_eof()        const noexcept {return file_.eof();}

    bool load_next(system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        std::ignore = detail::read_bytes_as<std::int32_t>(file_); // magic
        std::ignore = detail::read_bytes_as<std::int32_t>(file_); // version
        const auto len_title = detail::read_bytes_as<std::int32_t>(file_);
        std::vector<char> title(len_title);
        file_.read(title.data(), len_title);

        // TODO: skip unused blocks

        std::ignore = detail::read_bytes_as<std::int32_t>(file_);
        std::ignore = detail::read_bytes_as<std::int32_t>(file_);
        const auto boundary_size = detail::read_bytes_as<std::int32_t>(file_);
        std::ignore = detail::read_bytes_as<std::int32_t>(file_);
        std::ignore = detail::read_bytes_as<std::int32_t>(file_);
        std::ignore = detail::read_bytes_as<std::int32_t>(file_);
        std::ignore = detail::read_bytes_as<std::int32_t>(file_);
        this->position_block_size_ = detail::read_bytes_as<std::int32_t>(file_);
        this->velocity_block_size_ = detail::read_bytes_as<std::int32_t>(file_);
        this->force_block_size_    = detail::read_bytes_as<std::int32_t>(file_);
        this->number_of_particles_ = detail::read_bytes_as<std::int32_t>(file_);
        std::ignore = detail::read_bytes_as<std::int32_t>(file_); // current step
        std::ignore = detail::read_bytes_as<std::int32_t>(file_);
        std::ignore = detail::read_bytes_as<real_type>(file_);
        std::ignore = detail::read_bytes_as<real_type>(file_);

        this->has_unitcell_ = (boundary_size != 0);
        this->read_unitcell_if_needed(sys.boundary());

        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const real_type x = detail::read_bytes_as<real_type>(file_);
            const real_type y = detail::read_bytes_as<real_type>(file_);
            const real_type z = detail::read_bytes_as<real_type>(file_);
            sys.position(i) = sys.adjust_position(
                    math::make_coordinate<coordinate_type>(x, y, z));
        }
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const real_type x = detail::read_bytes_as<real_type>(file_);
            const real_type y = detail::read_bytes_as<real_type>(file_);
            const real_type z = detail::read_bytes_as<real_type>(file_);
            sys.velocity(i) = math::make_coordinate<coordinate_type>(x, y, z);
        }
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const real_type x = detail::read_bytes_as<real_type>(file_);
            const real_type y = detail::read_bytes_as<real_type>(file_);
            const real_type z = detail::read_bytes_as<real_type>(file_);
            sys.force(i) = math::make_coordinate<coordinate_type>(x, y, z);
        }

        return true;
    }

    std::string const& filename() const noexcept override {return filename_;}

  private:

    void read_unitcell_if_needed(
        UnlimitedBoundary<real_type, coordinate_type>&) noexcept
    {
        return; // No boundary exists. Do nothing.
    }
    void read_unitcell_if_needed(
        CuboidalPeriodicBoundary<real_type, coordinate_type>& bdry) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        if(!this->has_unitcell_)
        {
            MJOLNIR_LOG_WARN("TRR file lacks unit cell information");
            return;
        }

        // X axis vector
        const auto Xx = detail::read_bytes_as<real_type>(this->file_);
        const auto Xy = detail::read_bytes_as<real_type>(this->file_);
        const auto Xz = detail::read_bytes_as<real_type>(this->file_);

        // Y axis vector
        const auto Yx = detail::read_bytes_as<real_type>(this->file_);
        const auto Yy = detail::read_bytes_as<real_type>(this->file_);
        const auto Yz = detail::read_bytes_as<real_type>(this->file_);

        // Z axis vector
        const auto Zx = detail::read_bytes_as<real_type>(this->file_);
        const auto Zy = detail::read_bytes_as<real_type>(this->file_);
        const auto Zz = detail::read_bytes_as<real_type>(this->file_);

        if(Xy != real_type(0) || Xz != real_type(0) ||
           Yx != real_type(0) || Yz != real_type(0) ||
           Zx != real_type(0) || Zy != real_type(0))
        {
            MJOLNIR_LOG_WARN("The unit cell is not a rectangle");
            MJOLNIR_LOG_WARN("X axis vector = ", Xx, ", ", Xy, ", ", Xz);
            MJOLNIR_LOG_WARN("X axis vector = ", Yx, ", ", Yy, ", ", Yz);
            MJOLNIR_LOG_WARN("X axis vector = ", Zx, ", ", Zy, ", ", Zz);
        }

        // set width; XXX it requires sys.adjust_position later because the
        // original lower bound is unknown (only the width is recorded).
        const auto& lw = bdry.lower_bound();
        bdry.set_upper_bound(math::make_coordinate<coordinate_type>(
                    math::X(lw) + Xx, math::Y(lw) + Yy, math::Z(lw) + Zz));
        return ;
    }

  private:

    bool          has_unitcell_;
    std::string   filename_;
    std::size_t   number_of_frames_;
    std::size_t   number_of_particles_;
    std::size_t   position_block_size_;
    std::size_t   velocity_block_size_;
    std::size_t   force_block_size_;
    std::ifstream file_;
};

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class TRRLoader<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class TRRLoader<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class TRRLoader<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class TRRLoader<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif

#endif//MJOLNIR_CORE_TRR_LOADER_HPP
