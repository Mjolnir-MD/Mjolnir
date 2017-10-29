#ifndef JARNGREIPR_WRITE_AS_XYZ
#define JARNGREIPR_WRITE_AS_XYZ
#include <jarngreipr/model/CGChain.hpp>
#include <jarngreipr/io/XYZData.hpp>
#include <jarngreipr/io/XYZWrite.hpp>

namespace mjolnir
{

template<typename T>
void write_as_xyz(std::ostream& os, const CGChain<T>& chain)
{
    os << chain.size() << '\n';
    os << chain.front().kind() << '\n';
    for(const auto& bead : chain)
    {
        os << XYZLine<T>{bead->name, bead->position()} << '\n';
    }
    return;
}

} // mjolnir
#endif // JARNGREIPR_WRITE_AS_XYZ
