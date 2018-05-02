#ifndef JARNGREIPR_XYZ_FRAME_HPP
#define JARNGREIPR_XYZ_FRAME_HPP
#include <jarngreipr/xyz/XYZLine.hpp>
#include <vector>

namespace jarngreipr
{

template<typename realT, typename coordT>
struct XYZFrame
{
    typedef realT  real_type;
    typedef coordT coordinate_type;
    typedef XYZLine<realT, coordT> line_type;

    std::string            comment;
    std::vector<line_type> lines;
};

} // jarngreipr
#endif // JARNGREIPR_XYZ_FRAME_HPP
