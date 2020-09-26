#ifndef MJOLNIR_CORE_MSGPACK_SAVER_HPP
#define MJOLNIR_CORE_MSGPACK_SAVER_HPP
#include <mjolnir/util/macro.hpp>
#include <mjolnir/util/binary_io.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <mjolnir/core/System.hpp>
#include <fstream>
#include <iomanip>

namespace mjolnir
{

// Serialize System and RNG into the following MsgPack format.
// In the current implementation, the order should be preserved.
//
// system:
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
//
// rng:
// fixmap<1>{"internal_state": string}

template<typename traitsT>
class MsgPackSaver
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using system_type     = System<traits_type>;
    using attribute_type  = typename system_type::attribute_type;
    using rng_type        = RandomNumberGenerator<traits_type>;

  public:

    explicit MsgPackSaver(const std::string& filename_prefix)
      : prefix_(filename_prefix)
    {}
    ~MsgPackSaver() {}

    void save(const rng_type& rng)
    {
        // sometimes we need both x_rng1.msg and x_rng2.msg...
        return this->save(rng, "_rng");
    }
    void save(const rng_type& rng, const std::string& suffix)
    {
        constexpr std::uint8_t fixmap1_code = 0x81;
        this->buffer_.clear();
        this->buffer_.push_back(fixmap1_code);

        to_msgpack("internal_state");
        to_msgpack(rng.internal_state());

        const auto filename = prefix_ + suffix + std::string(".msg");
        std::ofstream ofs(filename);
        if(!ofs.good())
        {
            throw std::runtime_error("file open error: " + filename);
        }
        ofs.write(reinterpret_cast<const char*>(buffer_.data()),buffer_.size());
        ofs.close();
        return;
    }

    void save(const system_type& sys)
    {
        return this->save(sys, "_system");
    }
    void save(const system_type& sys, const std::string& suffix)
    {
        // 0b'1000'0011
        // 0x    8    3
        constexpr std::uint8_t fixmap3_code = 0x83;

        // ---------------------------------------------------------------------
        // clear buffer before writing to it.
        this->buffer_.clear();
        this->buffer_.push_back(fixmap3_code);

        // ---------------------------------------------------------------------
        // write boundary condition.

        to_msgpack("boundary");
        to_msgpack(sys.boundary());

        // ---------------------------------------------------------------------
        // append particles

        to_msgpack("particles");

        const auto num_particles = sys.size();
        if(num_particles < 16)
        {
            // 0b'1001'0000
            // 0x    9    0
            buffer_.push_back(std::uint8_t(num_particles) + std::uint8_t(0x90));
        }
        else if(num_particles < 65536)
        {
            const std::uint16_t size = num_particles;
            buffer_.push_back(std::uint8_t(0xdc));
            to_big_endian(size);
        }
        else
        {
            const std::uint32_t size = num_particles;
            buffer_.push_back(std::uint8_t(0xdd));
            to_big_endian(size);
        }

        constexpr std::uint8_t fixmap6_code = 0x86;
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            buffer_.push_back(fixmap6_code);
            to_msgpack("mass"         );
            to_msgpack(sys.mass(i)    );
            to_msgpack("position"     );
            to_msgpack(sys.position(i));
            to_msgpack("velocity"     );
            to_msgpack(sys.velocity(i));
            to_msgpack("force"        );
            to_msgpack(sys.force(i)   );
            to_msgpack("name"         );
            to_msgpack(sys.name(i)    );
            to_msgpack("group"        );
            to_msgpack(sys.group(i)   );
        }

        // ---------------------------------------------------------------------
        // write attributes list

        to_msgpack("attributes");
        to_msgpack(sys.attributes());

        // -------------------------------------------------------------------
        // overwrite .msg file by the current status

        const std::string filename = prefix_ + suffix + std::string(".msg");
        std::ofstream ofs(filename);
        if(!ofs.good())
        {
            throw std::runtime_error("file open error: " + filename);
        }
        ofs.write(reinterpret_cast<const char*>(buffer_.data()),buffer_.size());
        ofs.close();
        return;
    }

    std::string const& prefix() const noexcept {return prefix_;}

  private:

    void to_msgpack(const std::string& str)
    {
        constexpr std::uint8_t str8_code  = 0xd9;
        constexpr std::uint8_t str16_code = 0xda;
        constexpr std::uint8_t str32_code = 0xdb;

        // add a byte tag and length
        if(str.size() < 32)
        {
            // 0b'1010'0000
            // 0x    a    0
            buffer_.push_back(std::uint8_t(str.size()) + std::uint8_t(0xa0));
        }
        else if(str.size() <= 0xFF)
        {
            const std::uint8_t size = str.size();
            buffer_.push_back(str8_code);
            buffer_.push_back(size);
        }
        else if(str.size() <= 0xFFFF)
        {
            buffer_.push_back(str16_code);
            const std::uint16_t size = str.size();
            to_big_endian(size);
        }
        else
        {
            buffer_.push_back(str32_code);
            const std::uint32_t size = str.size();
            to_big_endian(size);
        }

        // write string body
        for(const char c : str)
        {
            buffer_.push_back(// write a char as a std::uint8_t
                    *reinterpret_cast<const std::uint8_t*>(std::addressof(c)));
        }
        return ;
    }

    void to_msgpack(const float& x)
    {
        constexpr std::uint8_t f32_code = 0xca;
        buffer_.push_back(f32_code);
        to_big_endian(x);
        return ;
    }
    void to_msgpack(const double& x)
    {
        constexpr std::uint8_t f64_code = 0xcb;
        buffer_.push_back(f64_code);
        to_big_endian(x);
        return ;
    }

    void to_msgpack(const coordinate_type& v)
    {
        // [float, float, float]
        // fixmap (3)
        // 0b'1001'0011
        // 0x    9    3
        constexpr std::uint8_t fixarray3_code = 0x93;
        buffer_.push_back(fixarray3_code);
        to_msgpack(math::X(v));
        to_msgpack(math::Y(v));
        to_msgpack(math::Z(v));
        return ;
    }

    void to_msgpack(const UnlimitedBoundary<real_type, coordinate_type>&)
    {
        // UnlimitedBoundary is represented as nil
        buffer_.push_back(std::uint8_t(0xc0));
        return ;
    }

    void to_msgpack(
        const CuboidalPeriodicBoundary<real_type, coordinate_type>& boundary)
    {
        constexpr std::uint8_t fixmap2_t = 0x82;
        buffer_.push_back(fixmap2_t);
        to_msgpack("lower");
        to_msgpack(boundary.lower_bound());
        to_msgpack("upper");
        to_msgpack(boundary.upper_bound());
        return ;
    }

    void to_msgpack(const attribute_type& attr)
    {
        constexpr std::uint8_t map16_code = 0xde;
        constexpr std::uint8_t map32_code = 0xdf;

        if(attr.size() < 16)
        {
            buffer_.push_back(std::uint8_t(attr.size()) + std::uint8_t(0x80));
        }
        else if(attr.size() < 65536)
        {
            buffer_.push_back(map16_code);
            const std::uint16_t size = attr.size();
            to_big_endian(size);
        }
        else
        {
            buffer_.push_back(map32_code);
            const std::uint32_t size = attr.size();
            to_big_endian(size);
        }

        for(const auto& keyval : attr)
        {
            // string -> real
            to_msgpack(keyval.first);
            to_msgpack(keyval.second);
        }
        return ;
    }

    template<typename T>
    void to_big_endian(const T& val)
    {
        const char* src = reinterpret_cast<const char*>(std::addressof(val));

#if defined(MJOLNIR_WITH_LITTLE_ENDIAN)
        // If the architecture uses little endian, we need to reverse bytes.
        std::reverse_copy(src, src + sizeof(T), std::back_inserter(buffer_));
#elif defined(MJOLNIR_WITH_BIG_ENDIAN)
        // If the architecture uses big endian, we don't need to do anything.
        std::copy(src, src + sizeof(T), std::back_inserter(buffer_));
#else
#  error "Unknown platform."
#endif
        return ;
    }

  private:

    std::string prefix_;
    std::vector<std::uint8_t> buffer_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class MsgPackSaver<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class MsgPackSaver<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class MsgPackSaver<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class MsgPackSaver<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif // MJOLNIR_CORE_MSGPACK_SAVER_HPP
