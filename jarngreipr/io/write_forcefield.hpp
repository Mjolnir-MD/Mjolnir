#ifndef JARNGREIPR_WRITE_FORCEFIELD_HPP
#define JARNGREIPR_WRITE_FORCEFIELD_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <jarngreipr/io/write_toml_value.hpp>

namespace jarngreipr
{

template<typename charT, typename traits>
void write_local_forcefield(
        std::basic_ostream<charT, traits>& os, const toml::Table& ff)
{
    os << "[[forcefields.local]]\n";
    os << "interaction = ";
    write_toml_value(os, mjolnir::toml_value_at(
                ff, "interaction", "write_local_forcefield"));
    os << '\n';

    os << "potential = ";
    write_toml_value(os, mjolnir::toml_value_at(
                ff, "potential", "write_local_forcefield"));
    os << '\n';

    os << "topology = ";
    write_toml_value(os, mjolnir::toml_value_at(
                ff, "topology", "write_local_forcefield"));
    os << '\n';

    const auto& params = mjolnir::toml_value_at(
        ff, "parameters", "write_local_forcefield").cast<toml::value_t::Array>();

    os << "parameters = [\n";
    for(const auto& p : params)
    {
        write_toml_value(os, p);
        os << ",\n";
    }
    os << "]\n";

    return ;
}

template<typename charT, typename traits>
void write_global_forcefield(
        std::basic_ostream<charT, traits>& os, const toml::Table& ff)
{
    os << "[[forcefields.global]]\n";
    for(const auto& kv : ff)
    {
        if(kv.first == "parameters") {continue;}
        os << kv.first << " = ";
        write_toml_value(os, kv.second);
        os << '\n';
    }

    const auto& params = mjolnir::toml_value_at(
        ff, "parameters", "write_global_forcefield").cast<toml::value_t::Array>();

    os << "parameters = [\n";
    for(const auto& p : params)
    {
        write_toml_value(os, p);
        os << ",\n";
    }
    os << "]\n";
    return ;
}

template<typename charT, typename traits>
void write_forcefield(
        std::basic_ostream<charT, traits>& os, const toml::Table& ff)
{
    os << "[[forcefields]]\n";

    if(ff.count("local") == 1)
    {
        for(const auto& local : ff["local"].cast<toml::value_t::Array>())
        {
            write_local_forcefield(os, local.cast<toml::value_t::Table>());
        }
    }
    if(ff.count("global") == 1)
    {
        for(const auto& local : ff["global"].cast<toml::value_t::Array>())
        {
            write_local_forcefield(os, local.cast<toml::value_t::Table>());
        }
    }
    return;
}

} // jarngreipr
#endif // JARNGREIPR_WRITE_FORCEFIELD_HPP
