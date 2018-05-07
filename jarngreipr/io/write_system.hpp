#ifndef JARNGREIPR_WRITE_SYSTEM_HPP
#define JARNGREIPR_WRITE_SYSTEM_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <iomanip>

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
    for(const auto& item : ps)
    {
        const auto& particle = item.cast<toml::value_t::Table>();
        const auto m = toml::get<double>(mjolnir::toml_value_at(
            particle, "mass", "write_system: particles"));
        const auto p = toml::get<std::array<double, 3>>(mjolnir::toml_value_at(
            particle, "position", "write_system: particles"));
        const auto v = toml::get<std::array<double, 3>>(mjolnir::toml_value_at(
            particle, "velocity", "write_system: particles"));

        os << "{mass=" << std::setw(7) << std::setprecision(2)
           << std::fixed << std::right << m;
        os << ",position=["
           << std::setw(9) << std::setprecision(4) << std::fixed << std::right
           << p[0] << ','
           << std::setw(9) << std::setprecision(4) << std::fixed << std::right
           << p[1] << ','
           << std::setw(9) << std::setprecision(4) << std::fixed << std::right
           << p[2];
        os << "],velocity=["
           << std::setw(9) << std::setprecision(4) << std::fixed << std::right
           << v[0] << ','
           << std::setw(9) << std::setprecision(4) << std::fixed << std::right
           << v[1] << ','
           << std::setw(9) << std::setprecision(4) << std::fixed << std::right
           << v[2];
        os << "]},\n";
    }
    os << "]\n";
    return ;
}

} // jarngreipr
#endif //JARNGREIPR_WRITE_PARTICLES_HPP
