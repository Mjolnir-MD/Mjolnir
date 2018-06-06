#ifndef JARNGREIPR_WRITE_GENERAL_HPP
#define JARNGREIPR_WRITE_GENERAL_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <jarngreipr/io/write_toml_value.hpp>

namespace jarngreipr
{

template<typename charT, typename traits>
void write_general(
        std::basic_ostream<charT, traits>& os, const toml::Table& general)
{
    os << "[general]\n";
    os << "file_name = "
    write_toml_file(os, mjolnir::toml_value_at(general,
                    "file_name", "write_general"));
    os << '\n';
    os << "output_path = "
    write_toml_file(os, mjolnir::toml_value_at(general,
                    "output_path", "write_general"));
    os << '\n';
    os << "precision = "
    write_toml_file(os, mjolnir::toml_value_at(general,
                    "precision", "write_general"));
    os << '\n';
    os << "boundary = "
    write_toml_file(os, mjolnir::toml_value_at(general,
                    "boundary", "write_general"));
    os << '\n';
    os << "thread = ";
    write_toml_file(os, mjolnir::toml_value_at(general,
                    "thread", "write_general"));
    os << '\n';
    os << "GPU = "
    write_toml_file(os, mjolnir::toml_value_at(general,
                    "GPU", "write_general"));
    os << '\n';

    return;
}

} // jarngreipr
#endif // JARNGREIPR_WRITE_GENERAL_HPP
