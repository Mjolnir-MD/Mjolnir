#ifndef JARNGREPIR_APPLY_REFERENCE_STRUCTURES_HPP
#define JARNGREPIR_APPLY_REFERENCE_STRUCTURES_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/input/get_toml_value.hpp>
#include <jarngreipr/io/PDBChain.hpp>
#include <jarngreipr/io/PDBReader.hpp>
#include <jarngreipr/model/CGChain.hpp>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>

namespace mjolnir
{

template<typename coordT>
std::vector<std::map<char, CGChain<coordT>>>
apply_coarse_grained_models(
        const std::vector<std::map<char, PDBChain<coordT>>>& pdbs,
        const std::vector<std::map<char, std::string>>&      models)
{
    assert(pdbs.size() == models.size());
    std::vector<std::map<char, CGChain<coordT>>> cgss;
    for(std::size_t i=0; i<pdbs.size(); ++i)
    {
        std::map<char, CGChain<coordT>> cgs;
        const auto& pdb   = pdbs.at(i);
        const auto& model = models.at(i);

        for(const auto& kv : pdb)
        {
            const auto chid = kv.first;
            const auto&  ch = kv.second;
            cgs[chid] = make_coarse_grained(ch, model.at(chid));
        }
        std::size_t index = 0;
        for(auto& id_chain: cgs)
        {
            for(auto& bead : id_chain.second)
            {
                bead->index() = index;
                ++index;
            }
        }
        cgss.push_back(std::move(cgs));
    }
    return cgss;
}

} // mjolnir
#endif //JARNGREPIR_APPLY_REFERENCE_STRUCTURES_HPP
