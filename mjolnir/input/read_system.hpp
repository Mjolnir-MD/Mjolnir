#ifndef MJOLNIR_INPUT_READ_SYSTEM_HPP
#define MJOLNIR_INPUT_READ_SYSTEM_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/format_nth.hpp>
#include <mjolnir/math/vector_util.hpp>
#include <mjolnir/input/read_table_from_file.hpp>
#include <mjolnir/input/read_path.hpp>
#include <mjolnir/input/utility.hpp>

namespace mjolnir
{
namespace detail
{

// reads boundary settings depending on the type information.
// the specific cases are implemented below.
template<typename boundaryT> struct read_boundary_impl;

template<typename realT, typename coordT>
struct read_boundary_impl<UnlimitedBoundary<realT, coordT>>
{
    static UnlimitedBoundary<realT, coordT>
    invoke(const toml::value& v)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("no boundary is set. unlimited");
        check_keys_available(v, {}); // no available keys
        return UnlimitedBoundary<realT, coordT>{};
    }
};

template<typename realT, typename coordT>
struct read_boundary_impl<CuboidalPeriodicBoundary<realT, coordT>>
{
    static CuboidalPeriodicBoundary<realT, coordT>
    invoke(const toml::value& boundary)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("shape of periodic boundary is cuboid");
        check_keys_available(boundary, {"upper"_s, "lower"_s});

        const auto upper = toml::find<coordT>(boundary, "upper");
        const auto lower = toml::find<coordT>(boundary, "lower");

        if(math::X(upper) <= math::X(lower) ||
           math::Y(upper) <= math::Y(lower) || math::Z(upper) <= math::Z(lower))
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_boundary: upper should be larger than lower",
                toml::find(boundary, "upper"), "upper boundary here",
                toml::find(boundary, "lower"), "lower boundary here"));
        }
        MJOLNIR_LOG_INFO("upper limit of the boundary = ", upper);
        MJOLNIR_LOG_INFO("lower limit of the boundary = ", lower);

        return CuboidalPeriodicBoundary<realT, coordT>(lower, upper);
    }
};

} // detail

template<typename traitsT>
typename traitsT::boundary_type
read_boundary(const toml::value& boundary)
{
    using boundary_t = typename traitsT::boundary_type;
    return detail::read_boundary_impl<boundary_t>::invoke(boundary);
}

// forward declaration.
template<typename traitsT>
System<traitsT> load_system_from_msgpack(const std::string& msg_file);

// It reads particles and other system-specific attributes (e.g. temperature)
template<typename traitsT>
System<traitsT>
read_system(const toml::value& root, const std::size_t N)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type       = typename traitsT::real_type;
    using coordinate_type = typename traitsT::coordinate_type;

    const auto& systems = toml::find(root, "systems");

    MJOLNIR_LOG_NOTICE("reading ", format_nth(N), " system ...");

    // check [[systems]] has `file_name = "checkpoint.msg"`.
    // If `file_name` has `.msg` file, load system status from the file.
    if(systems.at(N).contains("file_name"))
    {
        if(1 < systems.at(N).size())
        {
            MJOLNIR_LOG_WARN("When `file_name` is specified, "
                             "all other keys are ignored.");
            check_keys_available(systems.at(N), {"file_name"_s});
        }

        const auto& fname = toml::find<std::string>(systems, N, "file_name");
        if(file_extension_is(fname, ".msg"))
        {
            const auto& input_path = get_input_path();

            MJOLNIR_LOG_NOTICE("msgpack file specified. load system status from ",
                               input_path, fname);

            return load_system_from_msgpack<traitsT>(input_path + fname);
        }
    }

    const auto system = read_table_from_file(systems.at(N), "systems");

    check_keys_available(system, {"boundary_shape"_s, "attributes"_s, "particles"_s});

    const auto& boundary  = toml::find<toml::value>(system, "boundary_shape");
    const auto& particles = toml::find<toml::array>(system, "particles");

    System<traitsT> sys(particles.size(), read_boundary<traitsT>(boundary));

    for(const auto& attr : toml::find<toml::table>(system, "attributes"))
    {
        const real_type attribute = toml::get<real_type>(attr.second);
        sys.attribute(attr.first) = attribute;
        MJOLNIR_LOG_INFO("attribute.", attr.first, " = ", attribute);
    }

    // if there is no particle, return.
    if(particles.empty())
    {
        return sys;
    }

    {
        // if the first particle has velocity, mjolnir assumes that velocity is
        // already initialized (velocity generation not required).
        const auto& p = particles.front();
        sys.velocity_initialized() = (p.contains("vel") || p.contains("velocity"));
    }

    // A functor to find a value that corresponds to either of the key.
    // If both key exists, throw an error.
    const auto find_either = [](const toml::value& v, const std::string& key1,
                                const std::string& key2) -> const toml::value&
        {
            if(v.contains(key1) && v.contains(key2) != 0)
            {
                throw_exception<std::runtime_error>(toml::format_error("[error]"
                    " key duplicates.", v.at(key1), "here", v.at(key2),
                    "this conflicts with the above value definition"));
            }
            if(v.contains(key1)) {return v.at(key1);}
            if(v.contains(key2)) {return v.at(key2);}

            throw_exception<std::out_of_range>(toml::format_error(
                "both keys, \""_s + key1 + "\" and \""_s + key2 +
                "\", are not found."_s, v, "in this table"_s));
        };

    MJOLNIR_LOG_NOTICE(particles.size(), " particles are found.");
    for(std::size_t i=0; i<particles.size(); ++i)
    {
        using coord_type = coordinate_type;
        const auto& p = particles.at(i);

        check_keys_available(p, {"mass"_s, "m"_s, "position"_s, "pos"_s,
                "velocity"_s, "vel"_s, "name"_s, "group"_s});

        sys.mass(i)     = toml::get<real_type >(find_either(p, "m", "mass"));
        sys.position(i) = toml::get<coord_type>(find_either(p, "pos", "position"));
        sys.velocity(i) = math::make_coordinate<coordinate_type>(0, 0, 0);

        if(sys.velocity_initialized()) // it requires `velocity`.
        {
            sys.velocity(i) = toml::get<coord_type>(find_either(p, "vel", "velocity"));
        }
        else // if velocity_initialized is false, velocity should not be there.
        {
            if(p.contains("vel") || p.contains("velocity"))
            {
                throw_exception<std::runtime_error>(toml::format_error("[error]"
                    "read_system(): missing (or extraneous) velocity input",
                    particles.front(), "this does not have `velocity` field",
                    p, "but this has `velocity` field", {
                    "To specify the initial velocity, it is necessary to give "
                    "the value at all particles.",
                    "To generate the initial velocity with a random number, "
                    "all the `velocity` fields must be empty."
                    }));
            }
        }

        sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
        sys.name(i)  = toml::find_or<std::string>(p, "name",  "X"   );
        sys.group(i) = toml::find_or<std::string>(p, "group", "NONE");
        sys.rmass(i) = 1.0 / sys.mass(i);

        MJOLNIR_LOG_INFO("mass = ",        sys.mass(i),
                          ", position = ", sys.position(i),
                          ", velocity = ", sys.velocity(i),
                          ", force = ",    sys.force(i),
                          ", name = ",     sys.name(i),
                          ", group = ",    sys.group(i));
    }
    return sys;
}

// ===========================================================================
// Msgpack related functions

template<typename InputIterator>
std::uint8_t read_a_byte(InputIterator& src, InputIterator sentinel)
{
    assert(src != sentinel);
    const std::uint8_t byte = *src; ++src;
    return byte;
}

template<typename T, typename InputIterator>
T from_big_endian(InputIterator& src, InputIterator sentinel)
{
    using difference_type = typename std::iterator_traits<InputIterator>::difference_type;
    assert(static_cast<difference_type>(sizeof(T)) <= std::distance(src, sentinel));

    T val;
    char* dst = reinterpret_cast<char*>(std::addressof(val));

#if defined(MJOLNIR_WITH_LITTLE_ENDIAN)
    // If the architecture uses little endian, we need to reverse bytes.
    std::reverse_copy(src, src + sizeof(T), dst);
#elif defined(MJOLNIR_WITH_BIG_ENDIAN)
    // If the architecture uses big endian, we don't need to do anything.
    std::copy(src, src + sizeof(T), dst);
#else
#  error "Unknown platform."
#endif
    src += sizeof(T);
    return val;
}

// load a float.
template<typename T, typename InputIterator>
typename std::enable_if<std::is_same<T, float>::value, T>::type
from_msgpack(InputIterator& iter, InputIterator sentinel)
{
    using difference_type = typename std::iterator_traits<InputIterator>::difference_type;
    assert(static_cast<difference_type>(sizeof(T)) <= std::distance(iter, sentinel));

    const std::uint8_t byte = read_a_byte(iter, sentinel);
    if(byte != 0xca)
    {
        throw_exception<std::runtime_error>("expected tag 0xCA, but got ",
                                            std::uint32_t(byte));
    }
    return from_big_endian<float>(iter, sentinel);
}

// load a double.
template<typename T, typename InputIterator>
typename std::enable_if<std::is_same<T, double>::value, T>::type
from_msgpack(InputIterator& iter, InputIterator sentinel)
{
    using difference_type = typename std::iterator_traits<InputIterator>::difference_type;
    assert(static_cast<difference_type>(sizeof(T)) <= std::distance(iter, sentinel));

    const std::uint8_t byte = read_a_byte(iter, sentinel);
    if(byte != 0xcb)
    {
        throw_exception<std::runtime_error>("expected tag 0xCB, but got ",
                                            std::uint32_t(byte));
    }
    return from_big_endian<double>(iter, sentinel);
}

// load a string.
template<typename T, typename InputIterator>
typename std::enable_if<std::is_same<T, std::string>::value, T>::type
from_msgpack(InputIterator& iter, InputIterator sentinel)
{
    constexpr std::uint8_t str8_code  = 0xd9;
    constexpr std::uint8_t str16_code = 0xda;
    constexpr std::uint8_t str32_code = 0xdb;

    const std::uint8_t tag = read_a_byte(iter, sentinel);
    // fixstr code: 0b'101x'xxxx.
    // in the range [0b'1010'0000 = 0xa0, 0b'1011'1111 = 0xbf]
    if(0xa0 <= tag && tag <= 0xbf)
    {
        // 0x101x'xxxx
        // 0x0001'1111
        const std::uint8_t len = (tag & 0x1f);
        assert(len <= std::distance(iter, sentinel));

        std::string str(iter, iter + static_cast<std::size_t>(len));
        iter += str.size();
        return str;
    }
    else if(tag == str8_code)
    {
        const std::uint8_t len = read_a_byte(iter, sentinel);
        assert(len <= std::distance(iter, sentinel));

        std::string str(iter, iter + static_cast<std::size_t>(len));
        iter += str.size();
        return str;
    }
    else if(tag == str16_code)
    {
        const std::uint16_t len = from_big_endian<std::uint16_t>(iter, sentinel);
        assert(len <= std::distance(iter, sentinel));

        std::string str(iter, iter + static_cast<std::size_t>(len));
        iter += str.size();
        return str;
    }
    else if(tag == str32_code)
    {
        const std::uint32_t len = from_big_endian<std::uint32_t>(iter, sentinel);
        assert(len <= std::distance(iter, sentinel));

        std::string str(iter, iter + static_cast<std::size_t>(len));
        iter += str.size();
        return str;
    }
    else
    {
        throw_exception<std::runtime_error>("expected string, but got different"
                " type-tag: ", std::uint32_t(tag));
    }
}

template<typename InputIterator>
void check_msgpack_key(const std::string& expected,
                       InputIterator& iter, InputIterator sentinel)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    // It's called so many times, thus it does not make a log-scope.

    // check key before reading body
    const auto key = from_msgpack<std::string>(iter, sentinel);
    if(key != expected)
    {
        MJOLNIR_LOG_ERROR("invalid format in system .msg file. "
                          "expected \"", expected, "\", got ", key);
        throw_exception<std::runtime_error>("failed to load .msg file");
    }
    return;
}

template<typename realT, typename coordT, typename InputIterator>
void load_boundary_from_msgpack(UnlimitedBoundary<realT, coordT>&,
                                InputIterator& iter, InputIterator sentinel)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    // read nil. thats' all.
    const std::uint8_t tag = read_a_byte(iter, sentinel);
    if(tag != 0xc0)
    {
        MJOLNIR_LOG_ERROR("invalid format in system .msg file. "
                          "expected nil(", 0xc0, "), got ", tag);
        throw_exception<std::runtime_error>("failed to load .msg file");
    }
    return ;
}

template<typename realT, typename coordT, typename Iterator>
void load_boundary_from_msgpack(CuboidalPeriodicBoundary<realT, coordT>& bdry,
                                Iterator& iter, Iterator sentinel)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    // fixmap<2>{
    //   "lower": [float, float, float]
    //   "upper": [float, float, float]
    // }
    {
        constexpr std::uint8_t fixmap2_code = 0x82;
        const std::uint8_t tag = read_a_byte(iter, sentinel);
        if(tag != fixmap2_code)
        {
            MJOLNIR_LOG_ERROR("invalid format in system .msg file. "
                              "expected 0x82, got ", std::uint32_t(tag));
            throw_exception<std::runtime_error>("failed to load .msg file");
        }
    }

    coordT lower, upper;
    check_msgpack_key("lower", iter, sentinel);
    {
        constexpr std::uint8_t fixarray3_code = 0x93;
        const std::uint8_t tag = read_a_byte(iter, sentinel);
        if(tag != fixarray3_code)
        {
            MJOLNIR_LOG_ERROR("invalid format in system .msg file. "
                              "expected 0x93, got ", std::uint32_t(tag));
            throw_exception<std::runtime_error>("failed to load .msg file");
        }
    }
    math::X(lower) = from_msgpack<realT>(iter, sentinel);
    math::Y(lower) = from_msgpack<realT>(iter, sentinel);
    math::Z(lower) = from_msgpack<realT>(iter, sentinel);

    check_msgpack_key("upper", iter, sentinel);
    {
        constexpr std::uint8_t fixarray3_code = 0x93;
        const std::uint8_t tag = read_a_byte(iter, sentinel);
        if(tag != fixarray3_code)
        {
            MJOLNIR_LOG_ERROR("invalid format in system .msg file. "
                              "expected 0x93, got ", std::uint32_t(tag));
            throw_exception<std::runtime_error>("failed to load .msg file");
        }
    }
    math::X(upper) = from_msgpack<realT>(iter, sentinel);
    math::Y(upper) = from_msgpack<realT>(iter, sentinel);
    math::Z(upper) = from_msgpack<realT>(iter, sentinel);

    bdry.set_boundary(lower, upper);
    return ;
}

// The order should be preserved.
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
System<traitsT> load_system_from_msgpack(const std::string& msg_file)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type       = typename traitsT::real_type;
    using boundary_type   = typename traitsT::boundary_type;

    // -----------------------------------------------------------------------
    // open file and check its status

    std::ifstream ifs(msg_file, std::ios_base::binary);
    if(!ifs.good())
    {
        MJOLNIR_LOG_ERROR("file open error: ", msg_file);
        throw std::runtime_error("[error] file open error: " + msg_file);
    }

    // -----------------------------------------------------------------------
    // count file size

    const auto beg = ifs.tellg();
    ifs.seekg(0, std::ios::end);
    const auto end = ifs.tellg();
    const auto file_size = end - beg;
    ifs.seekg(beg);

    // -----------------------------------------------------------------------
    // load all the bytes into memory

    std::vector<std::uint8_t> content(static_cast<std::size_t>(file_size));
    ifs.read(reinterpret_cast<char*>(content.data()), file_size);
    MJOLNIR_LOG_INFO("file content read");

    auto iter = content.begin();
    const auto sentinel = content.end();
    {
        constexpr std::uint8_t fixmap3_code = 0x83;
        const std::uint8_t tag = read_a_byte(iter, sentinel);
        if(tag != fixmap3_code)
        {
            MJOLNIR_LOG_ERROR("invalid format in system .msg file.");
            MJOLNIR_LOG_ERROR("the root object should be a map with 3 elems.");
            throw_exception<std::runtime_error>("failed to load .msg file");
        }
    }

    // -----------------------------------------------------------------------
    // load boundary condition

    boundary_type boundary;
    check_msgpack_key("boundary", iter, sentinel);
    load_boundary_from_msgpack(boundary, iter, sentinel);
    MJOLNIR_LOG_INFO("boundary loaded");

    // -----------------------------------------------------------------------
    // load particles

    check_msgpack_key("particles", iter, sentinel);

    std::size_t num_particles = 0;
    {
        constexpr std::uint8_t array16_code = 0xdc;
        constexpr std::uint8_t array32_code = 0xdd;

        const std::uint8_t tag = read_a_byte(iter, sentinel);
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
                    from_big_endian<std::uint16_t>(iter, sentinel));
        }
        else if(tag == array32_code)
        {
            num_particles = static_cast<std::size_t>(
                    from_big_endian<std::uint32_t>(iter, sentinel));
        }
        else
        {
            MJOLNIR_LOG_ERROR("invalid format in .msg file.");
            MJOLNIR_LOG_ERROR("expected type tag ", array32_code,
                              "(0xDD), but got ", std::uint32_t(tag));
            throw_exception<std::runtime_error>("failed to load .msg file");
        }
    }
    MJOLNIR_LOG_INFO("there are ", num_particles, " paricles stored");
    System<traitsT> sys(num_particles, boundary);

    constexpr std::uint8_t fixmap6_code = 0x86;
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        {
            const std::uint8_t tag = read_a_byte(iter, sentinel);
            if(tag != fixmap6_code)
            {
                MJOLNIR_LOG_ERROR("invalid format in system .msg file.");
                MJOLNIR_LOG_ERROR("the particle object should be fixmap<6>.");
                MJOLNIR_LOG_ERROR("but the code is ", std::uint32_t(tag));
                throw_exception<std::runtime_error>("failed to load .msg file");
            }
        }
        // load mass
        check_msgpack_key("mass", iter, sentinel);
        sys.mass(i)  = from_msgpack<real_type>(iter, sentinel);
        sys.rmass(i) = real_type(1) / sys.mass(i);

        // load position
        check_msgpack_key("position", iter, sentinel);
        {
            // 0b1001'0011
            // 0x   9    3
            const std::uint8_t tag = read_a_byte(iter, sentinel);
            if(tag != 0x93)
            {
                throw_exception<std::runtime_error>("expected tag 0x93, but got ",
                                                    std::uint32_t(tag));
            }
            math::X(sys.position(i)) = from_msgpack<real_type>(iter, sentinel);
            math::Y(sys.position(i)) = from_msgpack<real_type>(iter, sentinel);
            math::Z(sys.position(i)) = from_msgpack<real_type>(iter, sentinel);
        }

        // load velocity
        check_msgpack_key("velocity", iter, sentinel);
        {
            const std::uint8_t tag = read_a_byte(iter, sentinel);
            if(tag != 0x93)
            {
                throw_exception<std::runtime_error>("expected tag 0x93, but got ",
                                                    std::uint32_t(tag));
            }
            math::X(sys.velocity(i)) = from_msgpack<real_type>(iter, sentinel);
            math::Y(sys.velocity(i)) = from_msgpack<real_type>(iter, sentinel);
            math::Z(sys.velocity(i)) = from_msgpack<real_type>(iter, sentinel);
        }
        // load force
        check_msgpack_key("force", iter, sentinel);
        {
            const std::uint8_t tag = read_a_byte(iter, sentinel);
            if(tag != 0x93)
            {
                throw_exception<std::runtime_error>("expected tag 0x93, but got ",
                                                    std::uint32_t(tag));
            }
            math::X(sys.force(i)) = from_msgpack<real_type>(iter, sentinel);
            math::Y(sys.force(i)) = from_msgpack<real_type>(iter, sentinel);
            math::Z(sys.force(i)) = from_msgpack<real_type>(iter, sentinel);
        }

        // load name
        check_msgpack_key("name", iter, sentinel);
        sys.name(i) = from_msgpack<std::string>(iter, sentinel);

        // load group
        check_msgpack_key("group", iter, sentinel);
        sys.group(i) = from_msgpack<std::string>(iter, sentinel);
    }

    // since velocity values are loaded from .msg file, we don't need to
    // re-initialize system.velocity by random numbers.
    sys.velocity_initialized() = true;

    // -----------------------------------------------------------------------
    // load attributes

    check_msgpack_key("attributes", iter, sentinel);

    std::size_t num_attributes = 0;
    {
        constexpr std::uint8_t map16_code = 0xde;
        constexpr std::uint8_t map32_code = 0xdf;

        const std::uint8_t tag = read_a_byte(iter, sentinel);
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
                    from_big_endian<std::uint16_t>(iter, sentinel));
        }
        else if(tag == map32_code)
        {
            num_attributes = static_cast<std::size_t>(
                    from_big_endian<std::uint32_t>(iter, sentinel));
        }
        else
        {
            MJOLNIR_LOG_ERROR("invalid format in .msg file.");
            MJOLNIR_LOG_ERROR("expected type tag ", map32_code,
                              "(0xDF), but got ", std::uint32_t(tag));
            throw_exception<std::runtime_error>("failed to load .msg file");
        }
    }
    for(std::size_t i=0; i<num_attributes; ++i)
    {
        const auto key = from_msgpack<std::string>(iter, sentinel);
        const auto val = from_msgpack<real_type>(iter, sentinel);
        sys.attribute(key) = val;
    }
    return sys;
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template System<SimulatorTraits<double, UnlimitedBoundary>       > read_system<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value& root, std::size_t N);
extern template System<SimulatorTraits<float,  UnlimitedBoundary>       > read_system<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value& root, std::size_t N);
extern template System<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_system<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value& root, std::size_t N);
extern template System<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_system<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value& root, std::size_t N);
#endif

}//mjolnir
#endif //MJOLNIR_READ_SYSTEM
