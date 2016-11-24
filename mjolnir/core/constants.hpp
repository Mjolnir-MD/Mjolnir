#ifndef MJOLNIR_CONSTANTS
#define MJOLNIR_CONSTANTS

namespace mjolnir
{

template<typename traitsT>
struct constants
{
    typedef typename traitsT::real_type real_type;
    constexpr static real_type tolerance = 1e-10; // XXX
};

}//mjolnir
#endif /* MJOLNIR_CONSTANTS */
