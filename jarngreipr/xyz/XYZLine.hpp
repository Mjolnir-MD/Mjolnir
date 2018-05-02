#ifndef JARNGREIPR_IO_XYZ_LINE_HPP
#define JARNGREIPR_IO_XYZ_LINE_HPP
#include <mjolnir/math/Vector.hpp>
#include <utility>
#include <ostream>
#include <istream>
#include <iomanip>
#include <sstream>
#include <string>

namespace jarngreipr
{

template<typename realT>
struct XYZLine
{
    typedef realT  real_type;
    typedef mjolnir::Vector<real_type, 3> coordinate_type;

    std::string     name;
    coordinate_type position;
};

template<typename realT>
inline XYZLine<realT>
make_xyz_line(std::string name, const mjolnir::Vector<realT, 3>& pos)
{
    return XYZLine<realT, coordT>{std::move(name), pos};
}

template<typename charT, typename traits, typename realT>
std::basic_ostream<charT, traits>& operator<<(
    std::basic_ostream<charT, traits>& os, const XYZLine<realT>& xyz)
{
    os << std::setw(6)  << std::left << xyz.name
       << std::fixed << std::showpoint << std::setprecision(5) << std::right
       << std::setw(10) << xyz.position[0] << ' '
       << std::setw(10) << xyz.position[1] << ' '
       << std::setw(10) << xyz.position[2];
    return os;
}

template<typename charT, typename traits, typename realT>
std::basic_istream<charT, traits>& operator>>(
    std::basic_istream<charT, traits>& is, XYZLine<realT>& xyz)
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
