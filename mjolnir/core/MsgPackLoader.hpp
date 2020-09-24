#ifndef MJOLNIR_CORE_MSGPACK_LOADER_HPP
#define MJOLNIR_CORE_MSGPACK_LOADER_HPP
#include <mjolnir/core/LoaderBase.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/util/logger.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace mjolnir
{

// fixmap<3> {
//     "boundary"     : fixmap<2>{"lower": [real, real, real],
//                                "upper": [real, real, real]}
//                   or nil,
//     "particles"    : array<N>[
//         fixmap<6>{
//             "mass"    : real,
//             "position": [real, real, real],
//             "velocity": [real, real, real],
//             "force"   : [real, real, real],
//             "name"    : string,
//             "group"   : string,
//         }, ...
//     ]
//     "attributres"  : map<N>{"temperature": real, ...},
// }
template<typename traitsT>
class MsgPackLoader final : public LoaderBase<traitsT>
{
  public:
    using base_type         = LoaderBase<traitsT>;
    using traits_type       = typename base_type::traits_type;
    using real_type         = typename base_type::real_type;
    using coordinate_type   = typename base_type::coordinate_type;
    using system_type       = typename base_type::system_type;

  public:

    explicit MsgPackLoader(const std::string& filename)
        : base_type(), filename_(filename), number_of_particles_(0),
          file_(filename_, std::ios::binary | std::ios::in)
    {
        if(!file_.good())
        {
            throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    " file open error: ", filename_);
        }
    }
    ~MsgPackLoader() override {}

    void initialize() override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // -------------------------------------------------------------------
        // check the first type tag

        {
            constexpr std::uint8_t fixmap3_code = 0x83;
            const auto tag = detail::read_bytes_as<std::uint8_t>(file_);
            if(tag != fixmap3_code)
            {
                MJOLNIR_LOG_ERROR("invalid format in system .msg file.");
                MJOLNIR_LOG_ERROR("the root object should be a map with 3 elems.");
                throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                        "failed to load ", filename_, " file");
            }
        }

        // -------------------------------------------------------------------
        // skip boundary
        {
            this->check_msgpack_key("boundary");
            typename system_type::boundary_type b;
            this->load_boundary(b);
        }
        MJOLNIR_LOG_INFO("boundary map skipped");

        // -------------------------------------------------------------------
        // check number of particles

        this->check_msgpack_key("particles");

        this->number_of_particles_ = 0;
        {
            constexpr std::uint8_t array16_code = 0xdc;
            constexpr std::uint8_t array32_code = 0xdd;

            const auto tag = detail::read_bytes_as<std::uint8_t>(file_);

            // fixarray size tag: 1001xxxx.
            // the possible values are ... [1001'0000 = 0x90, 1001'1111 = 0x9F]
            if(0x90 <= tag && tag <= 0x9F)
            {
                // take the last 4 bits
                this->number_of_particles_ = static_cast<std::size_t>(tag & 0x0F);
            }
            else if(tag == array16_code)
            {
                this->number_of_particles_ = static_cast<std::size_t>(
                        this->from_big_endian<std::uint16_t>());
            }
            else if(tag == array32_code)
            {
                this->number_of_particles_ = static_cast<std::size_t>(
                        this->from_big_endian<std::uint32_t>());
            }
            else
            {
                MJOLNIR_LOG_ERROR("invalid format in .msg file.");
                MJOLNIR_LOG_ERROR("expected type tag 0xDD, but got ",
                                  std::hex, std::uint32_t(tag));
                throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                        "failed to load ", filename_, " file");
            }
        }

        // -------------------------------------------------------------------
        // rewind and reset bitflags inside

        this->file_.clear();                      // clear bitflags
        this->file_.seekg(0, std::ios_base::beg); // rewind
        this->file_.peek();                       // update bitflags
        return;
    }

    std::size_t num_particles() const noexcept override {return number_of_particles_;}
    std::size_t num_frames()    const noexcept override {return 1;}
    bool        is_eof()        const noexcept override {return file_.eof();}

    bool load_next(system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // -------------------------------------------------------------------
        // check the first type tag
        {
            constexpr std::uint8_t fixmap3_code = 0x83;
            const auto tag = detail::read_bytes_as<std::uint8_t>(file_);
            if(tag != fixmap3_code)
            {
                MJOLNIR_LOG_ERROR("invalid format in system .msg file.");
                MJOLNIR_LOG_ERROR("the root object should be a map with 3 elems.");
                throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                        "failed to load ", filename_, " file");
            }
        }

        // -------------------------------------------------------------------
        // load boundary
        {
            this->check_msgpack_key("boundary");
            this->load_boundary(sys.boundary());
        }

        // -------------------------------------------------------------------
        // load number of particles

        check_msgpack_key("particles");

        this->number_of_particles_ = 0;
        {
            constexpr std::uint8_t array16_code = 0xdc;
            constexpr std::uint8_t array32_code = 0xdd;

            const auto tag = detail::read_bytes_as<std::uint8_t>(file_);

            // fixarray size tag: 1001xxxx.
            // the possible values are ... [1001'0000 = 0x90, 1001'1111 = 0x9F]
            if(0x90 <= tag && tag <= 0x9F)
            {
                // take the last 4 bits
                this->number_of_particles_ = static_cast<std::size_t>(tag & 0x0F);
            }
            else if(tag == array16_code)
            {
                this->number_of_particles_ = static_cast<std::size_t>(
                        this->from_big_endian<std::uint16_t>());
            }
            else if(tag == array32_code)
            {
                this->number_of_particles_ = static_cast<std::size_t>(
                        this->from_big_endian<std::uint32_t>());
            }
            else
            {
                MJOLNIR_LOG_ERROR("invalid format in .msg file.");
                MJOLNIR_LOG_ERROR("expected type tag 0xDD, but got ",
                                  std::hex, std::uint32_t(tag));
                throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                        "failed to load ", filename_, " file");
            }
        }
        if(sys.size() != this->number_of_particles_)
        {
            throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    " invalid number of particles. ", this->number_of_particles_,
                    " in the msg file, ", sys.size(), " in the system.");
        }

        // -------------------------------------------------------------------
        // load particles

        constexpr std::uint8_t fixmap6_code = 0x86;
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            {
                const auto tag = detail::read_bytes_as<std::uint8_t>(file_);
                if(tag != fixmap6_code)
                {
                    MJOLNIR_LOG_ERROR("invalid format in ", this->filename_,
                            ". the particle object should be fixmap<6> (0x86)."
                            " but the code is ", std::hex, std::uint32_t(tag));
                    throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                        "failed to load ", filename_, " file");
                }
            }

            // load mass
            this->check_msgpack_key("mass");
            this->from_msgpack(sys.mass(i));
            sys.rmass(i) = real_type(1) / sys.mass(i);

            // load position
            this->check_msgpack_key("position");
            {
                const auto tag = detail::read_bytes_as<std::uint8_t>(file_);
                if(tag != 0x93)
                {
                    throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                            "expected tag 0x93, but got ", std::hex, std::uint32_t(tag));
                }
                this->from_msgpack(math::X(sys.position(i)));
                this->from_msgpack(math::Y(sys.position(i)));
                this->from_msgpack(math::Z(sys.position(i)));
            }
            // load velocity
            this->check_msgpack_key("velocity");
            {
                const auto tag = detail::read_bytes_as<std::uint8_t>(file_);
                if(tag != 0x93)
                {
                    throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                            "expected tag 0x93, but got ", std::hex, std::uint32_t(tag));
                }
                this->from_msgpack(math::X(sys.velocity(i)));
                this->from_msgpack(math::Y(sys.velocity(i)));
                this->from_msgpack(math::Z(sys.velocity(i)));
            }
            // load force
            this->check_msgpack_key("force");
            {
                const auto tag = detail::read_bytes_as<std::uint8_t>(file_);
                if(tag != 0x93)
                {
                    throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                            "expected tag 0x93, but got ", std::hex, std::uint32_t(tag));
                }
                this->from_msgpack(math::X(sys.force(i)));
                this->from_msgpack(math::Y(sys.force(i)));
                this->from_msgpack(math::Z(sys.force(i)));
            }
            // load name
            this->check_msgpack_key("name");
            this->from_msgpack(sys.name(i));
            // load group
            this->check_msgpack_key("group");
            this->from_msgpack(sys.group(i));
        }

        // since velocity values are loaded from .msg file, we don't need to
        // re-initialize system.velocity by random numbers.
        sys.velocity_initialized() = true;

        // -----------------------------------------------------------------------
        // load attributes

        this->check_msgpack_key("attributes");

        std::size_t num_attributes = 0;
        {
            constexpr std::uint8_t map16_code = 0xde;
            constexpr std::uint8_t map32_code = 0xdf;
            const auto tag = detail::read_bytes_as<std::uint8_t>(file_);

            // fixmap size tag: 1000xxxx.
            // the possible values are ... [1000'0000 = 0x80, 1000'1111 = 0x8F]
            if(0x80 <= tag && tag <= 0x8F)
            {
                // take the last 4 bits
                num_attributes = static_cast<std::size_t>(tag & 0x0F);
            }
            else if(tag == map16_code)
            {
                num_attributes = static_cast<std::size_t>(
                        from_big_endian<std::uint16_t>());
            }
            else if(tag == map32_code)
            {
                num_attributes = static_cast<std::size_t>(
                        from_big_endian<std::uint32_t>());
            }
            else
            {
                MJOLNIR_LOG_ERROR("invalid format in .msg file.");
                MJOLNIR_LOG_ERROR("expected type tag ", map32_code,
                                  "(0xDF), but got ", std::uint32_t(tag));
                throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    "failed to load ", filename_, " file");
            }
        }
        for(std::size_t i=0; i<num_attributes; ++i)
        {
            std::string key;
            this->from_msgpack(key);

            real_type val;
            this->from_msgpack(val);

            sys.attribute(key) = val;
        }
        return true;
    }

    std::string const& filename() const noexcept override {return filename_;}

  private:

    void load_boundary(UnlimitedBoundary<real_type, coordinate_type>&)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // nil
        const auto tag = detail::read_bytes_as<std::uint8_t>(file_);
        if(tag != 0xc0)
        {
            MJOLNIR_LOG_ERROR("invalid format in system .msg file. expected "
                "nil tag (0xC0), got ", std::hex, std::uint32_t(tag));
            throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    "failed to load ", filename_, " file");
        }
        return ;
    }
    void load_boundary(CuboidalPeriodicBoundary<real_type, coordinate_type>& b)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // fixmap<2>{
        //   "lower": [float, float, float]
        //   "upper": [float, float, float]
        // }
        {
            constexpr std::uint8_t fixmap2_code = 0x82;
            const auto tag = detail::read_bytes_as<std::uint8_t>(file_);
            if(tag != fixmap2_code)
            {
                MJOLNIR_LOG_ERROR("invalid format in system .msg file. "
                    "expected 0x82, got ", std::hex, std::uint32_t(tag));
                throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    "failed to load ", filename_, " file");
            }
        }

        coordinate_type lower, upper;
        this->check_msgpack_key("lower");
        {
            constexpr std::uint8_t fixarray3_code = 0x93;
            const auto tag = detail::read_bytes_as<std::uint8_t>(file_);
            if(tag != fixarray3_code)
            {
                MJOLNIR_LOG_ERROR("invalid format in system .msg file. "
                                  "expected 0x93, got ", std::uint32_t(tag));
                throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    "failed to load ", filename_, " file");
            }
        }
        this->from_msgpack(math::X(lower));
        this->from_msgpack(math::Y(lower));
        this->from_msgpack(math::Z(lower));

        this->check_msgpack_key("upper");
        {
            constexpr std::uint8_t fixarray3_code = 0x93;
            const auto tag = detail::read_bytes_as<std::uint8_t>(file_);
            if(tag != fixarray3_code)
            {
                MJOLNIR_LOG_ERROR("invalid format in system .msg file. "
                                  "expected 0x93, got ", std::uint32_t(tag));
                throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    "failed to load ", filename_, " file");
            }
        }
        this->from_msgpack(math::X(upper));
        this->from_msgpack(math::Y(upper));
        this->from_msgpack(math::Z(upper));

        b.set_boundary(lower, upper);
        return ;
    }

    void check_msgpack_key(const std::string& expected)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        // It's called so many times, thus it does not make a log-scope.
        this->check_file_state();

        // check key before reading body
        std::string key;
        from_msgpack(key);
        if(key != expected)
        {
            MJOLNIR_LOG_ERROR("invalid format in system .msg file. "
                              "expected \"", expected, "\", got \"", key, "\"");
            throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    "failed to load ", filename_, " file");
        }
        MJOLNIR_LOG_INFO("key \"", expected, "\" found.");
        return;
    }

    void from_msgpack(std::string& str)
    {
        this->check_file_state();

        constexpr std::uint8_t str8_code  = 0xd9;
        constexpr std::uint8_t str16_code = 0xda;
        constexpr std::uint8_t str32_code = 0xdb;

        const auto tag = detail::read_bytes_as<std::uint8_t>(file_);
        // fixstr code: 0b'101x'xxxx.
        // in the range [0b'1010'0000 = 0xa0, 0b'1011'1111 = 0xbf]
        if(0xa0 <= tag && tag <= 0xbf)
        {
            // 0x101x'xxxx
            // 0x0001'1111
            const std::uint8_t len = (tag & 0x1f);
            std::copy_n(std::istream_iterator<char>(file_), len,
                        std::back_inserter(str));
            return;
        }
        else if(tag == str8_code)
        {
            const auto len = detail::read_bytes_as<std::uint8_t>(file_);;
            std::copy_n(std::istream_iterator<char>(file_), len,
                        std::back_inserter(str));
            return;
        }
        else if(tag == str16_code)
        {
            const auto len = from_big_endian<std::uint16_t>();
            std::copy_n(std::istream_iterator<char>(file_), len,
                        std::back_inserter(str));
            return;
        }
        else if(tag == str32_code)
        {
            const auto len = from_big_endian<std::uint32_t>();
            std::copy_n(std::istream_iterator<char>(file_), len,
                        std::back_inserter(str));
            return;
        }
        else
        {
            throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    "expected string (0xa0-0xbf, 0xd9, 0xda, 0xdb), but got "
                    "different type-tag: ", std::hex, std::uint32_t(tag));
        }
        return;
    }

    void from_msgpack(float& value)
    {
        this->check_file_state();

        const auto tag = detail::read_bytes_as<std::uint8_t>(file_);
        if(tag != 0xca)
        {
            throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                "expected float (type tag = 0xCA), but got ",
                std::hex, std::uint32_t(tag));
        }
        value = this->from_big_endian<float>();
        return;
    }
    void from_msgpack(double& value)
    {
        this->check_file_state();

        const auto tag = detail::read_bytes_as<std::uint8_t>(file_);
        if(tag != 0xcb)
        {
            throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                "expected double (type tag = 0xCB), but got ",
                std::hex, std::uint32_t(tag));
        }
        value = this->from_big_endian<double>();
        return;
    }

    template<typename T>
    T from_big_endian()
    {
        this->check_file_state();

        T val;
        char* dst = reinterpret_cast<char*>(std::addressof(val));

        std::array<char, sizeof(T)> buffer;
        this->file_.read(buffer.data(), buffer.size());

#if defined(MJOLNIR_WITH_LITTLE_ENDIAN)
        std::reverse_copy(buffer.begin(), buffer.end(), dst);
#elif defined(MJOLNIR_WITH_BIG_ENDIAN)
        std::copy(buffer.begin(), buffer.end(), dst);
#else
#  error "Unknown platform."
#endif
        return val;
    }

    void check_file_state()
    {
        if(file_.eof())
        {
            throw_exception<std::runtime_error>(
                "[error] mjolnir::MsgPackLoader: no more bytes to be read.");
        }
    }

  private:

    std::string   filename_;
    std::size_t   number_of_particles_;
    std::ifstream file_;
};

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class MsgPackLoader<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class MsgPackLoader<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class MsgPackLoader<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class MsgPackLoader<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif

#endif//MJOLNIR_CORE_MSGPACK_LOADER_HPP
