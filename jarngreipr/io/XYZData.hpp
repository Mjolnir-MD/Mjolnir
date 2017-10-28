#ifndef JARNGREIPR_IO_XYZ_DATA
#define JARNGREIPR_IO_XYZ_DATA
#include <utility>
#include <vector>
#include <string>

namespace mjolnir
{

template<typename coordT>
struct XYZLine
{
    typedef coordT coordinate_type;

    XYZLine()  = default;
    ~XYZLine() = default;
    XYZLine(const XYZLine&) = default;
    XYZLine(XYZLine&&)      = default;
    XYZLine& operator=(const XYZLine&) = default;
    XYZLine& operator=(XYZLine&&)      = default;

    explicit XYZLine(const std::string& str): name(str), position(0,0,0){}
    explicit XYZLine(std::string&& str): name(std::move(str)), position(0,0,0){}
    explicit XYZLine(const coordinate_type& crd): name(""), position(crd){}
    explicit XYZLine(coordinate_type&& crd): name(""), position(std::move(crd)){}

    XYZLine(const std::string& nm, const coordinate_type& crd)
        : name(nm), position(crd)
    {}
    XYZLine(std::string&& nm, const coordinate_type& crd)
        : name(std::move(nm)), position(crd)
    {}
    XYZLine(std::string&& nm, coordinate_type&& crd)
        : name(std::move(nm)), position(std::move(crd))
    {}

    std::string name;
    coordT      position;
};


template<typename coordT>
struct XYZFrame
{
    typedef coordT coordinate_type;
    typedef XYZLine<coordinate_type> line_type;

    XYZFrame()  = default;
    ~XYZFrame() = default;
    XYZFrame(const XYZFrame&) = default;
    XYZFrame(XYZFrame&&)      = default;
    XYZFrame& operator=(const XYZFrame&) = default;
    XYZFrame& operator=(XYZFrame&&)      = default;

    std::string            comment;
    std::vector<line_type> lines;
};

} // mjolnir
#endif //JARNGREIPR_IO_XYZ_DATA
