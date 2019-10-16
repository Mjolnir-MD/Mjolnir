#include <mjolnir/input/read_table_from_file.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template toml::basic_value<toml::discard_comments, std::unordered_map, std::vector>
read_table_from_file(const toml::basic_value<toml::discard_comments, std::unordered_map, std::vector>& root,
                     const std::string& name, const std::size_t N, const std::string& input_path);
} // mjolnir
