#ifndef JARNGREIPR_MAKE_COARSE_GRAINED
#define JARNGREIPR_MAKE_COARSE_GRAINED
#include <jarngreipr/model/CGChain.hpp>
#include <jarngreipr/io/PDBChain.hpp>

namespace mjolnir
{

template<typename coordT>
CGChain<coordT>
make_coarse_grained(const PDBChain<coordT>& pdb, const std::string& model)
{
    CGChain<coordT> cg;
    if(model == "CarbonAlpha")
    {
        for(const auto& residue : pdb)
        {
            cg.push_back(make_unique<CarbonAlpha<coordT>>(residue));
        }
    }
    else
    {
        throw std::runtime_error("jarngreipr::make_coarse_grained: "_str +
                "unknown model specified: "_str + model);
    }
}

} // mjolnir
#endif // JARNGREIPR_MAKE_COARSE_GRAINED
