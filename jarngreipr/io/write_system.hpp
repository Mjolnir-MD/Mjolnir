#ifndef JARNGREIPR_WRITE_SYSTEM_HPP
#define JARNGREIPR_WRITE_SYSTEM_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/util/get_toml_value.hpp>

namespace jarngreipr
{

template<typename charT, typename traits>
void write_system(std::basic_ostream<charT, traits>& os, const toml::Table& sys)
{
    os << "[[systems]]\n";
    os << "attributes = ";
    write_toml_value(os,
            mjolnir::toml_value_at(sys, "attributes", "write_system::sys"));
    os << '\n';

    os << "boundary = ";
    write_toml_value(os,
            mjolnir::toml_value_at(sys, "boundary", "write_system::sys"));
    os << '\n';

    const auto& ps = mjolnir::toml_value_at(
        sys, "particles", "write_system::sys").cast<toml::value_t::Array>();

    os << "particles = [\n";
    for(const auto& p : ps)
    {
        write_toml_value(os, p);
        os << ",\n";
    }
    os << "]\n";

    return ;
}

} // jarngreipr
#endif //JARNGREIPR_WRITE_PARTICLES_HPP
