#ifndef JARNGREIPR_IO_XYZ_WRITER
#define JARNGREIPR_IO_XYZ_WRITER
#include <jarngreipr/io/XYZData.hpp>
#include <stdexcept>
#include <ostream>
#include <fstream>
#include <sstream>

namespace mjolnir
{

template<typename coordT>
std::ostream& operator<<(std::ostream& ostrm, const XYZLine<coordT>& line)
{
    ostrm << line.name << std::fixed << std::showpoint << std::setprecision(5)
          << line.position[0] << ' ' << line.position[1] << ' '
          << line.position[2];
    return ostrm;
}

template<typename coordT>
std::ostream& operator<<(std::ostream& ostrm, const XYZFrame<coordT>& frame)
{
    ostrm << frame.lines.size() << '\n';
    ostrm << frame.comment      << '\n';
    for(const auto& line : frame.lines)
    {
        ostrm << line << '\n';
    }
    return ostrm;
}

template<typename coordT>
void
write_xyz_file(std::ostream& ostrm, const std::vector<XYZFrame<coordT>>& traj)
{
    for(const auto& frame : traj)
    {
        ostrm << frame;
    }
    return;
}

} // mjolnir
#endif// JARNGREIPR_IO_XYZ_WRITER
