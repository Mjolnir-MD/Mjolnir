#ifndef JARNGREIPR_XYZ_READER_HPP
#define JARNGREIPR_XYZ_READER_HPP
#include <jarngreipr/xyz/XYZLine.hpp>
#include <jarngreipr/xyz/XYZFrame.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <stdexcept>
#include <fstream>
#include <sstream>

namespace jarngreipr
{

// lazy XYZ file reader
template<typename realT>
class XYZReader
{
  public:
    typedef XYZLine<realT>  line_type;
    typedef XYZFrame<realT> frame_type;

  public:
    explicit XYZReader(const std::string& fname)
        : filename_(fname), ifstrm_(fname)
    {
        if(!ifstrm_.good())
        {
            throw std::runtime_error("jarngreipr::XYZReader: file open error: "
                    + filename_);
        }
    }
    ~XYZReader() = default;

    frame_type read_next_frame()
    {
        std::string line;
        std::getline(ifstrm_, line);
        std::size_t n;
        try
        {
            n = std::stoull(line);
        }
        catch(std::invalid_argument const& iva)
        {
            mjolnir::throw_exception<std::runtime_error>(
                "jarngreipr::XYZReader: in ", filename_,
                " invalid first line(not a number) appeared: ", line);
        }
        catch(std::out_of_range const& iva)
        {
            mjolnir::throw_exception<std::runtime_error>(
                "jarngreipr::XYZReader: in ", filename_,
                " invalid first line(not a number) appeared: ", line);
        }

        frame_type frame;
        std::getline(ifstrm_, frame.comment);
        for(std::size_t i=0; i<n; ++i)
        {
            std::getline(ifstrm_, line);
            std::istringstream iss(line);
            line_type xyz;
            iss >> xyz;
            frame.lines.push_back(xyz);
        }
        return frame;
    }

  private:
    std::string filename_;
    std::ifstream ifstrm_;
};

} // mjolnir
#endif //JARNGREIPR_IO_XYZ_READER
