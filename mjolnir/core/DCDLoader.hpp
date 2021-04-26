#ifndef MJOLNIR_CORE_DCD_LOADER_HPP
#define MJOLNIR_CORE_DCD_LOADER_HPP
#include <mjolnir/core/LoaderBase.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/util/logger.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace mjolnir
{

template<typename traitsT>
class DCDLoader final : public LoaderBase<traitsT>
{
  public:
    using base_type         = LoaderBase<traitsT>;
    using traits_type       = typename base_type::traits_type;
    using real_type         = typename base_type::real_type;
    using coordinate_type   = typename base_type::coordinate_type;
    using system_type       = typename base_type::system_type;

  public:

    explicit DCDLoader(const std::string& filename) : base_type(),
        has_unitcell_(false), filename_(filename),
        number_of_frames_(0), number_of_particles_(0),
        file_(filename_, std::ios::binary | std::ios::in)
    {
        if(!file_.good())
        {
            throw_exception<std::runtime_error>("[error] mjolnir::DCDLoader: "
                    "file open error: ", filename_);
        }
    }
    ~DCDLoader() override {}

    void initialize() override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // read header region and set `num_particles` and `num_frames`.

        // --------------------------------------------------------------------
        // the first block
        {
            const auto block_beg = detail::read_bytes_as<std::int32_t>(file_);

            std::ignore = detail::read_bytes_as<std::int32_t>(file_); // signature

            this->number_of_frames_ = detail::read_bytes_as<std::int32_t>(file_);

            std::ignore = detail::read_bytes_as<std::int32_t>(file_); // index_of_first
            std::ignore = detail::read_bytes_as<std::int32_t>(file_); // save_interval
            std::ignore = detail::read_bytes_as<std::int32_t>(file_); // total_step
            std::ignore = detail::read_bytes_as<std::int32_t>(file_); // total_chains

            detail::skip_bytes(file_, 4 * sizeof(std::int32_t)); // 4x int flags

            std::ignore = detail::read_bytes_as<float>(file_);
            this->has_unitcell_ = (detail::read_bytes_as<std::int32_t>(file_) == 1);

            detail::skip_bytes(file_, 8 * sizeof(std::int32_t)); // 8x int flags

            std::ignore = detail::read_bytes_as<std::int32_t>(file_);

            const auto block_end = detail::read_bytes_as<std::int32_t>(file_);

            MJOLNIR_LOG_NOTICE("There are ", this->number_of_frames_, " frames in ",
                               this->filename_);

            if(block_beg != block_end)
            {
                MJOLNIR_LOG_WARN("DCD file ", filename_, " seems to be broken.");
                MJOLNIR_LOG_WARN("Size of the first header block is inconsistent");
            }
        }

        // --------------------------------------------------------------------
        // the second block
        {
            const auto block_size_beg  = detail::read_bytes_as<std::int32_t>(file_);
            const auto number_of_lines = detail::read_bytes_as<std::int32_t>(file_);
            for(std::int32_t i=0; i<number_of_lines; ++i)
            {
                std::vector<char> comments(81, '\0');
                file_.read(comments.data(), 80);
                const std::string line(comments.begin(), comments.end());
                MJOLNIR_LOG_NOTICE("comment: ", line);
            }
            const auto block_size_end = detail::read_bytes_as<std::int32_t>(file_);

            if(block_size_beg != block_size_end)
            {
                MJOLNIR_LOG_WARN("DCD file ", filename_, " seems to be broken.");
                MJOLNIR_LOG_WARN("Size of the second header block is inconsistent");
            }
        }

        // --------------------------------------------------------------------
        // the third block
        {
            const auto block_size_beg  = detail::read_bytes_as<std::int32_t>(file_);
            this->number_of_particles_ = detail::read_bytes_as<std::int32_t>(file_);
            const auto block_size_end  = detail::read_bytes_as<std::int32_t>(file_);

            if(block_size_beg != block_size_end)
            {
                MJOLNIR_LOG_WARN("DCD file ", filename_, " seems to be broken.");
                MJOLNIR_LOG_WARN("Size of the second header block is inconsistent");
            }
        }
        MJOLNIR_LOG_NOTICE("There are ", this->number_of_particles_,
                           " particles in ", this->filename_);

        buffer_x_.resize(this->number_of_particles_);
        buffer_y_.resize(this->number_of_particles_);
        buffer_z_.resize(this->number_of_particles_);
        return;
    }

    std::size_t num_particles() const noexcept override {return number_of_particles_;}
    std::size_t num_frames()    const noexcept override {return number_of_frames_;}
    bool        is_eof()        const noexcept override {return file_.eof();}

    bool load_next(system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->read_unitcell_if_needed(sys.boundary());
        {
            const auto block_beg = detail::read_bytes_as<std::int32_t>(file_);

            assert(buffer_x_.size() == this->number_of_particles_);
            file_.read(reinterpret_cast<char*>(this->buffer_x_.data()),
                       this->number_of_particles_ * sizeof(float));

            const auto block_end = detail::read_bytes_as<std::int32_t>(file_);
            if(block_beg != block_end)
            {
                MJOLNIR_LOG_WARN("DCD file ", filename_, " seems to be broken.");
                MJOLNIR_LOG_WARN("Size of the X coordinate block is inconsistent."
                                 " Block header says there are ", block_beg,
                  " bytes, and the block footer says ", block_end, " bytes.");
                return false;
            }
        }
        {
            const auto block_beg = detail::read_bytes_as<std::int32_t>(file_);

            assert(buffer_y_.size() == this->number_of_particles_);
            file_.read(reinterpret_cast<char*>(this->buffer_y_.data()),
                       this->number_of_particles_ * sizeof(float));

            const auto block_end = detail::read_bytes_as<std::int32_t>(file_);
            if(block_beg != block_end)
            {
                MJOLNIR_LOG_WARN("DCD file ", filename_, " seems to be broken.");
                MJOLNIR_LOG_WARN("Size of the Y coordinate block is inconsistent."
                                 " Block header says there are ", block_beg,
                  " bytes, and the block footer says ", block_end, " bytes.");
                return false;
            }
        }
        {
            const auto block_beg = detail::read_bytes_as<std::int32_t>(file_);

            assert(buffer_z_.size() == this->number_of_particles_);
            file_.read(reinterpret_cast<char*>(this->buffer_z_.data()),
                       this->number_of_particles_ * sizeof(float));

            const auto block_end = detail::read_bytes_as<std::int32_t>(file_);
            if(block_beg != block_end)
            {
                MJOLNIR_LOG_WARN("DCD file ", filename_, " seems to be broken.");
                MJOLNIR_LOG_WARN("Size of the Z coordinate block is inconsistent."
                                 " Block header says there are ", block_beg,
                  " bytes, and the block footer says ", block_end, " bytes.");

                return false;
            }
        }

        for(std::size_t i=0; i<sys.size(); ++i)
        {
            sys.position(i) = sys.adjust_position(math::make_coordinate<
                coordinate_type>(buffer_x_[i], buffer_y_[i], buffer_z_[i]));
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
            MJOLNIR_LOG_WARN("dcd file ", this->filename_, " lacks unit cell "
                             "information.");
            return ;
        }

        const auto block_size_beg = detail::read_bytes_as<std::int32_t>(file_);
        const auto A              = detail::read_bytes_as<double>(file_);
        const auto gamma          = detail::read_bytes_as<double>(file_);
        const auto B              = detail::read_bytes_as<double>(file_);
        const auto beta           = detail::read_bytes_as<double>(file_);
        const auto alpha          = detail::read_bytes_as<double>(file_);
        const auto C              = detail::read_bytes_as<double>(file_);
        const auto block_size_end = detail::read_bytes_as<std::int32_t>(file_);

        if(block_size_beg != block_size_end)
        {
            MJOLNIR_LOG_WARN("DCD file ", filename_, " seems to be broken.");
            MJOLNIR_LOG_WARN("Size of the unitcell block is inconsistent");
        }

        constexpr double half_pi = math::constants<double>::half_pi();

        if((differs(alpha, 90.0) && differs(alpha, std::cos(half_pi))) ||
           (differs(beta,  90.0) && differs(beta,  std::cos(half_pi))) ||
           (differs(gamma, 90.0) && differs(gamma, std::cos(half_pi))))
        {
            MJOLNIR_LOG_WARN("The unit cell is not a rectangle");
            MJOLNIR_LOG_WARN("angle alpha = ", alpha);
            MJOLNIR_LOG_WARN("angle beta  = ", beta);
            MJOLNIR_LOG_WARN("angle gamma = ", gamma);
        }

        // set width; XXX it requires sys.adjust_position later because the
        // original lower bound is unknown (only the width is recorded).
        const auto& lw = bdry.lower_bound();
        bdry.set_upper_bound(math::make_coordinate<coordinate_type>(
                    math::X(lw) + A, math::Y(lw) + B, math::Z(lw) + C));
        return ;
    }

    static bool differs(const double lhs, const double rhs) noexcept
    {
        constexpr double rel_tol = math::rel_tolerance<double>();
        constexpr double abs_tol = math::abs_tolerance<double>();

        return abs_tol                     < std::abs(lhs - rhs) ||
               rel_tol * 0.5 * (lhs + rhs) < std::abs(lhs - rhs);
    }

  private:

    bool               has_unitcell_;
    std::string        filename_;
    std::size_t        number_of_frames_;
    std::size_t        number_of_particles_;
    std::vector<float> buffer_x_;
    std::vector<float> buffer_y_;
    std::vector<float> buffer_z_;
    std::ifstream      file_;
};

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class DCDLoader<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class DCDLoader<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class DCDLoader<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class DCDLoader<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif

#endif//MJOLNIR_CORE_DCD_LOADER_HPP
