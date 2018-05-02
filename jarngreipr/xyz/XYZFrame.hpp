#ifndef JARNGREIPR_XYZ_FRAME_HPP
#define JARNGREIPR_XYZ_FRAME_HPP
#include <jarngreipr/xyz/XYZLine.hpp>
#include <vector>

namespace jarngreipr
{

template<typename realT>
struct XYZFrame
{
    typedef realT  real_type;
    typedef XYZLine<realT> line_type;

    std::string            comment;
    std::vector<line_type> lines;
};

} // jarngreipr
#endif // JARNGREIPR_XYZ_FRAME_HPP
