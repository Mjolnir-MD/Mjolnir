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
    if(params.empty()) {return;}

    /* ============ prepair some values for formating parameters ============ */
    std::size_t max_index = 0;
    for(const auto& p : params)
    {
        const auto& ptable = p.cast<toml::value_t::Table>();
        if(ptable.count("indices"))
        {
            const auto idxs = toml::get<std::vector<std::size_t>>(
                mjolnir::toml_value_at(ptable, "indices",
                    "jarngreipr: generated local forcefield"));

            max_index = std::max(max_index,
                    *std::max_element(idxs.begin(), idxs.end()));
        }
    }
    const std::size_t idx_width = std::to_string(max_index).size();

    std::vector<std::string> keys;
    for(const auto& kv : params.front().cast<toml::value_t::Table>())
    {
        const auto& key = kv.first;
        if(key != "indices")
        {
            keys.push_back(key);
        }
    }
    // and sort keys in alphabetical order to make the order uniform
    std::sort(keys.begin(), keys.end());
    /* ================================ end! ================================ */

    os << "parameters = [\n";
    for(const auto& p : params)
    {
        const auto& ptable = p.cast<toml::value_t::Table>();
        os << '{';
        // write indices first
        if(ptable.count("indices"))
        {
            const auto idxs = toml::get<std::vector<std::size_t>>(
                mjolnir::toml_value_at(ptable, "indices",
                    "jarngreipr: generated local forcefield"));

            os << "indices = [";
            for(const auto& idx : idxs)
            {
                os << std::setw(idx_width) << idx << ',';
            }
            os << "], ";
        }

        for(const auto& key : keys)
        {
            const auto& val = mjolnir::toml_value_at(ptable, key,
                    "jarngreipr: generated local forcefield");
            os << key << " = ";
            write_toml_value(os, val);
            os << ", ";
        }
        os << "},\n";
    }
    os << "]\n";

    return ;
}

template<typename charT, typename traits>
void write_global_forcefield(
        std::basic_ostream<charT, traits>& os, const toml::Table& ff)
{
    os << "[[forcefields.global]]\n";

    {
        std::vector<std::string> keys;
        for(const auto& kv : ff)
        {
            if(kv.first != "parameters") {keys.push_back(kv.first);}
        }
        std::sort(keys.begin(), keys.end());

        for(const auto& key : keys)
        {
            os << key << " = ";
            write_toml_value(os, mjolnir::toml_value_at(ff, key,
                    "jarngreipr: generated global forcefield"));
            os << '\n';
        }
    }

    {
        const auto& params = mjolnir::toml_value_at(
            ff, "parameters", "write_global_forcefield"
            ).cast<toml::value_t::Array>();
        if(params.empty()) {return;}

        /* ========== prepair some values for formating parameters ========== */
        std::size_t max_index = 0;
        for(const auto& p : params)
        {
            const auto& ptable = p.cast<toml::value_t::Table>();
            if(ptable.count("index"))
            {
                const auto idx = toml::get<std::size_t>(
                    mjolnir::toml_value_at(ptable, "index",
                        "jarngreipr: generated local forcefield"));
                max_index = std::max(max_index, idx);
            }
        }
        const std::size_t idx_width = std::to_string(max_index).size();

        std::vector<std::string> keys;
        for(const auto& kv : params.front().cast<toml::value_t::Table>())
        {
            const auto& key = kv.first;
            if(key != "index")
            {
                keys.push_back(key);
            }
        }
        // and sort keys in alphabetical order to make the order uniform
        std::sort(keys.begin(), keys.end());
        /* ============================== end! ============================== */

        os << "parameters = [\n";
        for(const auto& p : params)
        {
            const auto& ptable = p.cast<toml::value_t::Table>();
            os << '{';
            // write index first
            if(ptable.count("index"))
            {
                const auto idx = toml::get<std::size_t>(mjolnir::toml_value_at(
                    ptable, "index", "jarngreipr: generated local forcefield"));
                os << "index = "<< std::setw(idx_width) << idx << ", ";
            }

            for(const auto& key : keys)
            {
                const auto& val = mjolnir::toml_value_at(ptable, key,
                        "jarngreipr: generated local forcefield");
                os << key << " = ";
                write_toml_value(os, val);
                os << ", ";
            }
            os << "},\n";
        }
        os << "]\n";
    }
    return ;
}

template<typename charT, typename traits>
void write_forcefield(
        std::basic_ostream<charT, traits>& os, const toml::Table& ff)
{
    os << "[[forcefields]]\n";

    if(ff.count("local") == 1)
    {
        for(const auto& local : ff.at("local").cast<toml::value_t::Array>())
        {
            write_local_forcefield(os, local.cast<toml::value_t::Table>());
        }
    }
    if(ff.count("global") == 1)
    {
        for(const auto& global : ff.at("global").cast<toml::value_t::Array>())
        {
            write_global_forcefield(os, global.cast<toml::value_t::Table>());
        }
    }
    return;
}

} // jarngreipr
#endif // JARNGREIPR_WRITE_FORCEFIELD_HPP
