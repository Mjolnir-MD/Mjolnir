#ifndef JARNGREIPR_IO_XYZ_LINE_HPP
#define JARNGREIPR_IO_XYZ_LINE_HPP
#include <utility>
#include <ostream>
#include <istream>
#include <iomanip>
#include <sstream>
#include <string>

namespace jarngreipr
{

template<typename realT, typename coordT>
struct XYZLine
{
    typedef realT  real_type;
    typedef coordT coordinate_type;

    std::string     name;
    coordinate_type position;
};

template<typename realT, typename coordT>
XYZLine<realT, coordT> make_xyz_line(std::string name, const coordT& pos)
{
    return XYZLine<realT, coordT>{std::move(name), pos};
}

template<typename charT, typename traits, typename realT, typename coordT>
std::basic_ostream<charT, traits>& operator<<(
    std::basic_ostream<charT, traits>& os, const XYZLine<realT, coordT>& xyz)
{
    os << std::setw(6)  << std::left << xyz.name
       << std::fixed << std::showpoint << std::setprecision(5) << std::right
       << std::setw(10) << xyz.position[0] << ' '
       << std::setw(10) << xyz.position[1] << ' '
       << std::setw(10) << xyz.position[2];
    return os;
}

template<typename charT, typename traits, typename realT, typename coordT>
std::basic_istream<charT, traits>& operator>>(
    std::basic_istream<charT, traits>& is, XYZLine<realT, coordT>& xyz)
{
    std::string line;
    std::getline(is, line);
    std::istringstream iss(line);

    iss >> xyz.name;
    iss >> xyz.position[0];
    iss >> xyz.position[1];
    iss >> xyz.position[2];
    return is;
}

} // jarngreipr
#endif //JARNGREIPR_IO_XYZ_DATA
