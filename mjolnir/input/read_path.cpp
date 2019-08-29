#include <mjolnir/input/read_path.hpp>

#ifndef MJOLNIR_SEPARATE_BUILD
#error "MJOLNIR_SEPARATE_BUILD flag is required"
#endif

namespace mjolnir
{
template std::string read_input_path <toml::discard_comments, std::unordered_map, std::vector>(const toml::basic_value<toml::discard_comments, std::unordered_map, std::vector>&);
template std::string read_output_path<toml::discard_comments, std::unordered_map, std::vector>(const toml::basic_value<toml::discard_comments, std::unordered_map, std::vector>&);
} // mjolnir
