#ifndef MJOLNIR_CORE_MSGPACK_OBSERVER_HPP
#define MJOLNIR_CORE_MSGPACK_OBSERVER_HPP
#include <mjolnir/util/macro.hpp>
#include <mjolnir/core/ObserverBase.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/Unit.hpp>
#include <fstream>
#include <iomanip>

namespace mjolnir
{

// Serialize System into MsgPack format that is equivalent to the following JSON
// In the current implementation, the order should be preserved.
// fixmap<3> {
//     "boundary"     : {"lower": [float, float, float],
//                       "upper": [float, float, float]} or nil,
//     "particles"    : [fixmap<7>{
//          "mass"    : float,
//          "rmass"   : float,
//          "position": [float, float, float],
//          "velocity": [float, float, float],
//          "force"   : [float, float, float],
//          "name"    : string,
//          "group"   : string,
//          }, ...
//     ]
//     "attributres"  : {"temperature": float, ...},
// }

template<typename traitsT>
class MsgPackObserver final : public ObserverBase<traitsT>
{
  public:
    using base_type         = ObserverBase<traitsT>;
    using traits_type       = typename base_type::traits_type;
    using real_type         = typename base_type::real_type;
    using coordinate_type   = typename base_type::coordinate_type;
    using system_type       = typename base_type::system_type;
    using forcefield_type   = typename base_type::forcefield_type;

    // system attribute container type. map-like class.
    using attribute_type    = typename system_type::attribute_type;

  public:

    explicit MsgPackObserver(const std::string& filename_prefix)
      : base_type(), prefix_(filename_prefix),
        filename_(filename_prefix + std::string(".msg"))
    {}
    ~MsgPackObserver() override {}

    void initialize(const std::size_t,  const real_type,
                    const system_type&, const forcefield_type&) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_NOTICE("checkpoint file is ", filename_);

        // check the specified file can be opened.
        // Here, it does not clear the content.
        std::ofstream ofs(filename_, std::ios_base::app);
        if(!ofs.good())
        {
            MJOLNIR_LOG_ERROR("file open error: ", filename_);
            throw std::runtime_error("file open error");
        }
        ofs.close();
        return;
    }

    void update(const std::size_t,  const real_type,
                const system_type&, const forcefield_type&) override
    {
        return; // do nothing.
    }

    void output(const std::size_t, const real_type,
                const system_type& sys, const forcefield_type&) override
    {
        // 0b'1000'0011
        // 0x    8    3
        constexpr std::uint8_t fixmap3_code = 0x83;

        // ---------------------------------------------------------------------
        // clear buffer before writing to it.
        this->buffer_.clear();
        this->buffer_.push_back(fixmap3_code);

        auto buffer_iterator = std::back_inserter(buffer_);

        // ---------------------------------------------------------------------
        // write boundary condition.

        to_msgpack("boundary",     buffer_iterator);
        to_msgpack(sys.boundary(), buffer_iterator);

        // ---------------------------------------------------------------------
        // append particles

        to_msgpack("particles", buffer_iterator);

        const auto num_particles = sys.size();
        if(num_particles < 16)
        {
            // 0b'1001'0000
            // 0x    9    0
            std::uint8_t fixarray_code = num_particles;
            fixarray_code |= std::uint8_t(0x90);
            buffer_.push_back(fixarray_code);
        }
        else if(num_particles < 65536)
        {
            const std::uint16_t size = num_particles;
            buffer_.push_back(std::uint8_t(0xdc));
            to_big_endian(size, buffer_iterator);
        }
        else
        {
            const std::uint32_t size = num_particles;
            buffer_.push_back(std::uint8_t(0xdd));
            to_big_endian(size, buffer_iterator);
        }

        constexpr std::uint8_t fixmap6_code = 0x86;
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            buffer_.push_back(fixmap6_code);
            to_msgpack("mass",          buffer_iterator);
            to_msgpack(sys.mass(i),     buffer_iterator);
            to_msgpack("position",      buffer_iterator);
            to_msgpack(sys.position(i), buffer_iterator);
            to_msgpack("velocity",      buffer_iterator);
            to_msgpack(sys.velocity(i), buffer_iterator);
            to_msgpack("force",         buffer_iterator);
            to_msgpack(sys.force(i),    buffer_iterator);
            to_msgpack("name",          buffer_iterator);
            to_msgpack(sys.name(i),     buffer_iterator);
            to_msgpack("group",         buffer_iterator);
            to_msgpack(sys.group(i),    buffer_iterator);
        }

        // ---------------------------------------------------------------------
        // write attributes list

        to_msgpack("attributes",     buffer_iterator);
        to_msgpack(sys.attributes(), buffer_iterator);

        // -------------------------------------------------------------------
        // overwrite .msg file by the current status

        std::ofstream ofs(filename_);
        if(!ofs.good())
        {
            throw std::runtime_error("file open error: " + filename_);
        }
        ofs.write(reinterpret_cast<const char*>(buffer_.data()),buffer_.size());
        ofs.close();
        return;
    }

    void finalize(const std::size_t step, const real_type t,
                  const system_type& sys, const forcefield_type& ff) override
    {
        this->output(step, t, sys, ff);
        return;
    }

    std::string const& prefix() const noexcept override {return prefix_;}

  private:

    template<typename OutputIterator>
    void to_msgpack(const std::string& str, OutputIterator& out)
    {
        constexpr std::uint8_t str8_code  = 0xd9;
        constexpr std::uint8_t str16_code = 0xda;
        constexpr std::uint8_t str32_code = 0xdb;

        // add a byte tag and length
        if(str.size() < 32)
        {
            // 0b'1010'0000
            // 0x    a    0
            std::uint8_t fixstr_code = str.size();
            fixstr_code |= std::uint8_t(0xa0);
            *out = fixstr_code; ++out;
        }
        else if(str.size() <= 0xFF)
        {
            const std::uint8_t size = str.size();
            *out = str8_code; ++out;
            *out = size;      ++out;
        }
        else if(str.size() <= 0xFFFF)
        {
            const std::uint16_t size = str.size();
            *out = str16_code; ++out;
            to_big_endian(size, out);
        }
        else
        {
            const std::uint32_t size = str.size();
            *out = str32_code; ++out;
            to_big_endian(size, out);
        }

        // write string body
        for(const std::uint8_t c : str)
        {
            *out = c; ++out;
        }
        return ;
    }

    template<typename OutputIterator>
    void to_msgpack(const float& x, OutputIterator& out)
    {
        constexpr std::uint8_t f32_code = 0xca;
        *out = f32_code; ++out;
        to_big_endian(x, out);
        return ;
    }
    template<typename OutputIterator>
    void to_msgpack(const double& x, OutputIterator& out)
    {
        constexpr std::uint8_t f64_code = 0xcb;
        *out = f64_code; ++out;
        to_big_endian(x, out);
        return ;
    }

    template<typename OutputIterator>
    void to_msgpack(const coordinate_type& v, OutputIterator& out)
    {
        // [float, float, float]
        // fixmap (3)
        // 0b'1001'0011
        // 0x    9    3
        constexpr std::uint8_t fixarray3_code = 0x93;
        *out = fixarray3_code; ++out;
        to_msgpack(math::X(v), out);
        to_msgpack(math::Y(v), out);
        to_msgpack(math::Z(v), out);
        return ;
    }

    template<typename OutputIterator>
    void to_msgpack(const UnlimitedBoundary<real_type, coordinate_type>&,
                    OutputIterator& out)
    {
        // UnlimitedBoundary is represented as nil
        *out = std::uint8_t(0xc0); ++out;
        return ;
    }

    template<typename OutputIterator>
    void to_msgpack(
        const CuboidalPeriodicBoundary<real_type, coordinate_type>& boundary,
        OutputIterator& out)
    {
        constexpr std::uint8_t fixmap2_t = 0x82;
        *out = fixmap2_t; ++out;
        to_msgpack("lower", out);
        to_msgpack(boundary.lower_bound(), out);
        to_msgpack("upper", out);
        to_msgpack(boundary.upper_bound(), out);
        return ;
    }

    template<typename OutputIterator>
    void to_msgpack(const attribute_type& attr, OutputIterator& out)
    {
        constexpr std::uint8_t map16_code = 0xde;
        constexpr std::uint8_t map32_code = 0xdf;

        if(attr.size() < 16)
        {
            std::uint8_t fixmap_code = attr.size();
            fixmap_code |= std::uint8_t(0x80);
            *out = fixmap_code; ++out;
        }
        else if(attr.size() < 65536)
        {
            const std::uint16_t size = attr.size();
            *out = map16_code; ++out;
            to_big_endian(size, out);
        }
        else
        {
            const std::uint32_t size = attr.size();
            *out = map32_code; ++out;
            to_big_endian(size, out);
        }

        for(const auto& keyval : attr)
        {
            // string -> real
            to_msgpack(keyval.first,  out);
            to_msgpack(keyval.second, out);
        }
        return ;
    }

    template<typename T, typename OutputIterator>
    void to_big_endian(const T& val, OutputIterator& dst)
    {
        const char* src = reinterpret_cast<const char*>(std::addressof(val));

#if defined(MJOLNIR_WITH_LITTLE_ENDIAN)
        // If the architecture uses little endian, we need to reverse bytes.
        std::reverse_copy(src, src + sizeof(T), dst);
#elif defined(MJOLNIR_WITH_BIG_ENDIAN)
        // If the architecture uses big endian, we don't need to do anything.
        std::copy(src, src + sizeof(T), dst);
#else
#  error "Unknown platform."
#endif
        return ;
    }

  private:

    std::string prefix_;
    std::string filename_;
    std::vector<std::uint8_t> buffer_;
};
#ifdef MJOLNIR_SEPARATE_BUILD
extern template class MsgPackObserver<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class MsgPackObserver<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class MsgPackObserver<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class MsgPackObserver<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif // MJOLNIR_CORE_MSGPACK_OBSERVER_HPP
