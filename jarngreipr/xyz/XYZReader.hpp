#ifndef JARNGREIPR_IO_XYZ_READER
#define JARNGREIPR_IO_XYZ_READER
#include <jarngreipr/io/XYZData.hpp>
#include <jarngreipr/util/string.hpp>
#include <stdexcept>
#include <istream>
#include <fstream>
#include <sstream>

namespace mjolnir
{

template<typename coordT>
XYZLine<coordT> read_xyz_line(const std::string& line)
{
    XYZLine<coordT> l;
    std::istringstream iss(line);
    iss >> l.name >> l.position[0] >> l.position[1] >> l.position[2];
    return l;
}

template<typename coordT>
XYZFrame<coordT> read_xyz_frame(std::istream& istrm)
{
    XYZFrame<coordT> f;
    std::string line;

    std::getline(istrm, line);
    std::size_t N;
    try
    {
        N = std::stoull(line);
    }
    catch(...)
    {
        throw std::runtime_error("mjolnir::io::read_xyz_frame: "_str +
                "found invalid line: "_str + line);
    }
    f.lines.reserve(N);

    std::getline(istrm, f.comment);
    while(!istrm.eof() || N != 0)
    {
        std::getline(istrm, line);
        f.push_back(read_xyz_line(line));
        istrm.peek();
        --N;
    }
    if(N != 0)
    {
        throw std::runtime_error("mjolnir::io::read_xyz_frame: "
                "eof while reading a frame");
    }
    return f;
}

template<typename coordT>
std::vector<XYZFrame<coordT>> read_xyz_file(std::istream& istrm)
{
    std::vector<XYZFrame<coordT>> fs;
    try
    {
        while(!istrm.eof())
        {
            fs.push_back(read_xyz_frame(istrm));
            istrm.peek();
        }
    }
    catch(const std::runtime_error& rte)
    {
        std::cerr << "mjolnir::io::read_xyz_file: Error in reading frames; ";
        std::cerr << "read only first " << fs.size() << "frames\n";
        std::cerr << "What: " << rte.what() << std::endl;
    }
    return fs;
}

// lazy XYZ file reader
template<typename coordT>
struct XYZReader
{
    typedef XYZLine<coordT>  line_type;
    typedef XYZFrame<coordT> frame_type;

    explicit XYZReader(const std::string& fname): ifstrm_(fname){}
    ~XYZReader() = default;

    frame_type read_next_frame()
    {
        return read_xyz_frame(this->ifstrm_);
    }

  private:
    std::ifstream ifstrm_;
};

} // mjolnir
#endif //JARNGREIPR_IO_XYZ_READER
