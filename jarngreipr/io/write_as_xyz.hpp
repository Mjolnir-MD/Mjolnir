#ifndef JARNGREIPR_WRITE_AS_XYZ
#define JARNGREIPR_WRITE_AS_XYZ
#include <jarngreipr/model/CGChain.hpp>
#include <jarngreipr/io/XYZData.hpp>
#include <jarngreipr/io/XYZWriter.hpp>

namespace mjolnir
{

template<typename T>
void write_as_xyz(std::ostream& os, const CGChain<T>& chain)
{
    os << chain.size() << '\n';
    os << chain.front()->kind() << '\n';
    for(std::size_t i=0; i<chain.size(); ++i)
    {
        const auto& bead = chain.at(i);
        const XYZLine<T> xyz(bead->name(), bead->position());
        os << xyz << '\n';
    }
    return;
}

} // mjolnir
#endif // JARNGREIPR_WRITE_AS_XYZ
