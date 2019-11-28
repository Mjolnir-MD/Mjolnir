#ifndef MJOLNIR_CORE_XYZ_LOADER_HPP
#define MJOLNIR_CORE_XYZ_LOADER_HPP
#include <mjolnir/core/LoaderBase.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/util/logger.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace mjolnir
{

template<typename traitsT>
class XYZLoader final : public LoaderBase<traitsT>
{
  public:
    using base_type         = LoaderBase<traitsT>;
    using traits_type       = typename base_type::traits_type;
    using real_type         = typename base_type::real_type;
    using coordinate_type   = typename base_type::coordinate_type;
    using system_type       = typename base_type::system_type;

  public:

    explicit XYZLoader(const std::string& filename)
        : base_type(), filename_(filename), file_(filename_), line_number_(0),
          number_of_frames_(0), number_of_particles_(0)
    {
        if(!file_.good())
        {
            throw_exception<std::runtime_error>("[error] mjolnir::XYZLoader: "
                    "file open error: ", filename_);
        }
    }
    ~XYZLoader() override {}

    void initialize() override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        std::string buf;
        std::getline(this->file_, buf);
        const std::size_t N = std::stoull(buf);
        this->number_of_particles_ = N;

        MJOLNIR_LOG_INFO("the first line is \"", buf, "\", N = ", N);

        std::size_t number_of_lines = 1;
        while(std::getline(this->file_, buf))
        {
            ++number_of_lines;
        }
        if(number_of_lines % (N+2) != 0)
        {
            MJOLNIR_LOG_WARN("XYZ file ", this->filename_, " may contain an "
                             "incomplete frame or the number of particles "
                             "changes. This may cause incorrect progress bar "
                             "and incorrect number of steps to run.");
        }
        this->number_of_frames_ = number_of_lines / (N+2);
        MJOLNIR_LOG_NOTICE("There are ", this->number_of_frames_, " frames in ",
                           this->filename_);

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
        this->file_.peek(); // to check the file reaches to EOF or not.
        if(this->file_.eof())
        {
            return false;
        }

        std::string line;
        std::getline(file_, line);
        ++line_number_;
        std::size_t N_tmp;
        try
        {
            N_tmp = std::stoull(line);
        }
        catch(const std::exception&)
        {
            std::cout << " " << line_number_ << " | " << line << std::endl;
            throw;
        }
        const std::size_t N = N_tmp;

        if(sys.size() != N)
        {
            throw_exception<std::runtime_error>("[error] mjolnir::XYZLoader: "
                "The number of particles in the system differs from the xyz file",
                filename_);
        }

        // read a comment line. just skip this.
        std::getline(file_, line);
        ++line_number_;

        for(std::size_t i=0; i<N; ++i)
        {
            std::getline(file_, line);
            ++line_number_;

            std::istringstream iss(line);
            std::string name;
            real_type x, y, z;
            iss >> name >> x >> y >> z;
            if(iss.fail())
            {
                throw_exception<std::runtime_error>("[error] mjolnir::XYZLoader"
                    ": failed to load a snapshot\n", line_number_, " | ", line);
            }
            sys.position(i) = math::make_coordinate<coordinate_type>(x, y, z);
        }
        this->file_.peek();
        return true;
    }

    std::string const& filename() const noexcept override {return filename_;}

  private:

    std::string   filename_;
    std::ifstream file_;
    std::size_t   line_number_;
    std::size_t   number_of_frames_;
    std::size_t   number_of_particles_;
};

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class XYZLoader<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class XYZLoader<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class XYZLoader<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class XYZLoader<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif

#endif//MJOLNIR_CORE_XYZ_LOADER_HPP
