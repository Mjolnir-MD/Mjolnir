#ifndef JARNGREIPR_XYZ_WRITER_HPP
#define JARNGREIPR_XYZ_WRITER_HPP
#include <jarngreipr/xyz/XYZLine.hpp>
#include <jarngreipr/xyz/XYZFrame.hpp>
#include <stdexcept>
#include <ostream>
#include <fstream>
#include <sstream>

namespace jarngreipr
{

template<typename realT>
class XYZWriter
{
  public:
    typedef XYZLine<realT>  line_type;
    typedef XYZFrame<realT> frame_type;

  public:
    explicit XYZWriter(const std::string& fname)
        : filename_(fname), ofstrm_(fname)
    {
        if(!ofstrm_.good())
        {
            throw std::runtime_error("jarngreipr::XYZWriter: file open error: "
                    + filename_);
        }
    }
    ~XYZWriter() = default;

    void write_frame(const frame_type& frame)
    {
        ofstrm_ << frame.lines.size() << '\n';
        ofstrm_ << frame.comment      << '\n';
        for(const auto& line : frame.lines)
        {
            ofstrm_ << line << '\n';
        }
        return;
    }

  private:
    std::string filename_;
    std::ofstream ofstrm_;
};

} // mjolnir
#endif// JARNGREIPR_IO_XYZ_WRITER
