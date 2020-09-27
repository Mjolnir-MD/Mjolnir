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
class MsgPackLoader
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using boundary_type   = typename traits_type::boundary_type;
    using system_type     = System<traits_type>;

  public:

    MsgPackLoader() {}
    ~MsgPackLoader() {}

    system_type load_system(const std::string& filename)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        std::ifstream file(filename);
        if(!file.good())
        {
            MJOLNIR_LOG_ERROR("file open error: ", filename);
            throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    "failed to load ", filename, " file");

        }

        // -------------------------------------------------------------------
        // check the first type tag
        {
            constexpr std::uint8_t fixmap3_code = 0x83;
            const auto tag = detail::read_bytes_as<std::uint8_t>(file);
            if(tag != fixmap3_code)
            {
                MJOLNIR_LOG_ERROR("invalid format in system .msg file.");
                MJOLNIR_LOG_ERROR("the root object should be a map with 3 elems.");
                throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                        "failed to load ", filename, " file");
            }
        }

        // -------------------------------------------------------------------
        // load boundary

        boundary_type boundary{};
        this->check_msgpack_key(file, filename, "boundary");
        this->load_boundary(file, filename, boundary);

        // -------------------------------------------------------------------
        // load number of particles

        this->check_msgpack_key(file, filename, "particles");
        std::size_t num_particles = 0;
        {
            constexpr std::uint8_t array16_code = 0xdc;
            constexpr std::uint8_t array32_code = 0xdd;

            const auto tag = detail::read_bytes_as<std::uint8_t>(file);

            // fixarray size tag: 1001xxxx.
            // the possible values are ... [1001'0000 = 0x90, 1001'1111 = 0x9F]
            if(0x90 <= tag && tag <= 0x9F)
            {
                // take the last 4 bits
                num_particles = static_cast<std::size_t>(tag & 0x0F);
            }
            else if(tag == array16_code)
            {
                num_particles = static_cast<std::size_t>(
                        this->from_big_endian<std::uint16_t>(file, filename));
            }
            else if(tag == array32_code)
            {
                num_particles = static_cast<std::size_t>(
                        this->from_big_endian<std::uint32_t>(file, filename));
            }
            else
            {
                MJOLNIR_LOG_ERROR("invalid format in .msg file.");
                MJOLNIR_LOG_ERROR("expected type tag 0xDD, but got ",
                                  std::hex, std::uint32_t(tag));
                throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                        "failed to load ", filename, " file");
            }
        }
        system_type sys(num_particles, boundary);

        // -------------------------------------------------------------------
        // load particles

        constexpr std::uint8_t fixmap6_code = 0x86;
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            {
                const auto tag = detail::read_bytes_as<std::uint8_t>(file);
                if(tag != fixmap6_code)
                {
                    MJOLNIR_LOG_ERROR("invalid format in ", filename,
                            ". the particle object should be fixmap<6> (0x86)."
                            " but the code is ", std::hex, std::uint32_t(tag));
                    throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                        "failed to load ", filename, " file");
                }
            }

            // load mass
            this->check_msgpack_key(file, filename, "mass");
            sys.mass(i)  = this->from_msgpack<real_type>(file, filename);
            sys.rmass(i) = real_type(1) / sys.mass(i);

            // load position
            this->check_msgpack_key(file, filename, "position");
            {
                const auto tag = detail::read_bytes_as<std::uint8_t>(file);
                if(tag != 0x93)
                {
                    throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                            "expected tag 0x93, but got ", std::hex, std::uint32_t(tag));
                }
                math::X(sys.position(i)) = this->from_msgpack<real_type>(file, filename);
                math::Y(sys.position(i)) = this->from_msgpack<real_type>(file, filename);
                math::Z(sys.position(i)) = this->from_msgpack<real_type>(file, filename);
            }
            // load velocity
            this->check_msgpack_key(file, filename, "velocity");
            {
                const auto tag = detail::read_bytes_as<std::uint8_t>(file);
                if(tag != 0x93)
                {
                    throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                            "expected tag 0x93, but got ", std::hex, std::uint32_t(tag));
                }
                math::X(sys.velocity(i)) = this->from_msgpack<real_type>(file, filename);
                math::Y(sys.velocity(i)) = this->from_msgpack<real_type>(file, filename);
                math::Z(sys.velocity(i)) = this->from_msgpack<real_type>(file, filename);
            }
            // load force
            this->check_msgpack_key(file, filename, "force");
            {
                const auto tag = detail::read_bytes_as<std::uint8_t>(file);
                if(tag != 0x93)
                {
                    throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                            "expected tag 0x93, but got ", std::hex, std::uint32_t(tag));
                }
                math::X(sys.force(i)) = this->from_msgpack<real_type>(file, filename);
                math::Y(sys.force(i)) = this->from_msgpack<real_type>(file, filename);
                math::Z(sys.force(i)) = this->from_msgpack<real_type>(file, filename);
            }
            // load name
            this->check_msgpack_key(file, filename, "name");
            sys.name(i) = this->from_msgpack<std::string>(file, filename);
            // load group
            this->check_msgpack_key(file, filename, "group");
            sys.group(i) = this->from_msgpack<std::string>(file, filename);
        }

        // since velocity values are loaded from .msg file, we don't need to
        // re-initialize system.velocity by random numbers.
        sys.velocity_initialized() = true;

        // -----------------------------------------------------------------------
        // load attributes

        this->check_msgpack_key(file, filename, "attributes");

        std::size_t num_attributes = 0;
        {
            constexpr std::uint8_t map16_code = 0xde;
            constexpr std::uint8_t map32_code = 0xdf;
            const auto tag = detail::read_bytes_as<std::uint8_t>(file);

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
                    this->from_big_endian<std::uint16_t>(file, filename));
            }
            else if(tag == map32_code)
            {
                num_attributes = static_cast<std::size_t>(
                    this->from_big_endian<std::uint32_t>(file, filename));
            }
            else
            {
                MJOLNIR_LOG_ERROR("invalid format in .msg file.");
                MJOLNIR_LOG_ERROR("expected type tag ", map32_code,
                                  "(0xDF), but got ", std::uint32_t(tag));
                throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    "failed to load ", filename, " file");
            }
        }
        for(std::size_t i=0; i<num_attributes; ++i)
        {
            const auto key = this->from_msgpack<std::string>(file, filename);
            const auto val = this->from_msgpack<real_type>(file, filename);
            sys.attribute(key) = val;
        }
        return sys;
    }

  private:

    void load_boundary(std::ifstream& file, const std::string filename,
                       UnlimitedBoundary<real_type, coordinate_type>&)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // nil
        const auto tag = detail::read_bytes_as<std::uint8_t>(file);
        if(tag != 0xc0)
        {
            MJOLNIR_LOG_ERROR("invalid format in system .msg file. expected "
                "nil tag (0xC0), got ", std::hex, std::uint32_t(tag));
            throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    "failed to load ", filename, " file");
        }
        return ;
    }
    void load_boundary(std::ifstream& file, const std::string filename,
                       CuboidalPeriodicBoundary<real_type, coordinate_type>& b)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // fixmap<2>{
        //   "lower": [float, float, float]
        //   "upper": [float, float, float]
        // }
        {
            constexpr std::uint8_t fixmap2_code = 0x82;
            const auto tag = detail::read_bytes_as<std::uint8_t>(file);
            if(tag != fixmap2_code)
            {
                MJOLNIR_LOG_ERROR("invalid format in system .msg file. "
                    "expected 0x82, got ", std::hex, std::uint32_t(tag));
                throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    "failed to load ", filename, " file");
            }
        }

        coordinate_type lower, upper;
        this->check_msgpack_key(file, filename, "lower");
        {
            constexpr std::uint8_t fixarray3_code = 0x93;
            const auto tag = detail::read_bytes_as<std::uint8_t>(file);
            if(tag != fixarray3_code)
            {
                MJOLNIR_LOG_ERROR("invalid format in system .msg file. "
                                  "expected 0x93, got ", std::uint32_t(tag));
                throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    "failed to load ", filename, " file");
            }
        }
        math::X(lower) = this->from_msgpack<real_type>(file, filename);
        math::Y(lower) = this->from_msgpack<real_type>(file, filename);
        math::Z(lower) = this->from_msgpack<real_type>(file, filename);

        this->check_msgpack_key(file, filename, "upper");
        {
            constexpr std::uint8_t fixarray3_code = 0x93;
            const auto tag = detail::read_bytes_as<std::uint8_t>(file);
            if(tag != fixarray3_code)
            {
                MJOLNIR_LOG_ERROR("invalid format in system .msg file. "
                                  "expected 0x93, got ", std::uint32_t(tag));
                throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    "failed to load ", filename, " file");
            }
        }
        math::X(upper) = this->from_msgpack<real_type>(file, filename);
        math::Y(upper) = this->from_msgpack<real_type>(file, filename);
        math::Z(upper) = this->from_msgpack<real_type>(file, filename);

        b.set_boundary(lower, upper);
        return ;
    }

    void check_msgpack_key(std::ifstream& file, const std::string& filename,
                           const std::string& expected)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        // It's called so many times, thus it does not make a log-scope.
        this->check_file_state(file, filename);

        // check key before reading body
        const auto key = this->from_msgpack<std::string>(file, filename);
        if(key != expected)
        {
            MJOLNIR_LOG_ERROR("invalid format in system .msg file. "
                              "expected \"", expected, "\", got \"", key, "\"");
            throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    "failed to load key ", expected, " from ", filename, ".");
        }
        MJOLNIR_LOG_INFO("key \"", expected, "\" found.");
        return;
    }

    template<typename T>
    typename std::enable_if<std::is_same<T, std::string>::value, std::string>::type
    from_msgpack(std::ifstream& file, const std::string& filename)
    {
        this->check_file_state(file, filename);

        constexpr std::uint8_t str8_code  = 0xd9;
        constexpr std::uint8_t str16_code = 0xda;
        constexpr std::uint8_t str32_code = 0xdb;

        std::string str;
        const auto tag = detail::read_bytes_as<std::uint8_t>(file);
        // fixstr code: 0b'101x'xxxx.
        // in the range [0b'1010'0000 = 0xa0, 0b'1011'1111 = 0xbf]
        if(0xa0 <= tag && tag <= 0xbf)
        {
            // 0x101x'xxxx
            // 0x0001'1111
            const std::uint8_t len = (tag & 0x1f);
            std::copy_n(std::istream_iterator<char>(file), len,
                        std::back_inserter(str));
        }
        else if(tag == str8_code)
        {
            const auto len = detail::read_bytes_as<std::uint8_t>(file);;
            std::copy_n(std::istream_iterator<char>(file), len,
                        std::back_inserter(str));
        }
        else if(tag == str16_code)
        {
            const auto len = from_big_endian<std::uint16_t>(file, filename);
            std::copy_n(std::istream_iterator<char>(file), len,
                        std::back_inserter(str));
        }
        else if(tag == str32_code)
        {
            const auto len = from_big_endian<std::uint32_t>(file, filename);
            std::copy_n(std::istream_iterator<char>(file), len,
                        std::back_inserter(str));
        }
        else
        {
            throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    "expected string (0xa0-0xbf, 0xd9, 0xda, 0xdb), but got "
                    "different type-tag: ", std::hex, std::uint32_t(tag));
        }
        return str;
    }

    template<typename T>
    typename std::enable_if<std::is_same<T, float>::value, float>::type
    from_msgpack(std::ifstream& file, const std::string& filename)
    {
        this->check_file_state(file, filename);

        const auto tag = detail::read_bytes_as<std::uint8_t>(file);
        if(tag != 0xca)
        {
            throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                "expected float (type tag = 0xCA), but got ",
                std::hex, std::uint32_t(tag));
        }
        return this->from_big_endian<float>(file, filename);
    }
    template<typename T>
    typename std::enable_if<std::is_same<T, double>::value, double>::type
    from_msgpack(std::ifstream& file, const std::string& filename)
    {
        this->check_file_state(file, filename);

        const auto tag = detail::read_bytes_as<std::uint8_t>(file);
        if(tag != 0xcb)
        {
            throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                "expected double (type tag = 0xCB), but got ",
                std::hex, std::uint32_t(tag));
        }
        return this->from_big_endian<double>(file, filename);
    }

    template<typename T>
    T from_big_endian(std::ifstream& file, const std::string& filename)
    {
        this->check_file_state(file, filename);

        T val;
        char* dst = reinterpret_cast<char*>(std::addressof(val));

        std::array<char, sizeof(T)> buffer;
        file.read(buffer.data(), buffer.size());

#if defined(MJOLNIR_WITH_LITTLE_ENDIAN)
        std::reverse_copy(buffer.begin(), buffer.end(), dst);
#elif defined(MJOLNIR_WITH_BIG_ENDIAN)
        std::copy(buffer.begin(), buffer.end(), dst);
#else
#  error "Unknown platform."
#endif
        return val;
    }

    void check_file_state(std::ifstream& file, const std::string& filename)
    {
        if(file.eof())
        {
            throw_exception<std::runtime_error>("[error] mjolnir::MsgPackLoader:"
                    " no more bytes to be read in ", filename, ".");
        }
        return;
    }
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
