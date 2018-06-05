#ifndef JARNGREIPR_MODEL_CG_MODEL_GENERATOR_HPP
#define JARNGREIPR_MODEL_CG_MODEL_GENERATOR_HPP
#include <jarngreipr/model/CGChain.hpp>
#include <jarngreipr/pdb/PDBChain.hpp>
#include <mjolnir/util/throw_exception.hpp>

namespace jarngreipr
{

template<typename realT>
class CGModelGeneratorBase
{
  public:
    typedef realT real_type;
    typedef CGChain<real_type>  cg_chain_type;
    typedef PDBChain<real_type> pdb_chain_type;

  public:

    virtual ~CGModelGeneratorBase() = default;

    virtual cg_chain_type
    generate(const pdb_chain_type& pdb, const std::size_t offset) const = 0;
};

// CGModelGenerator manages several subclasses of CGModelGeneratorBase.
// this is not a superclass. this is a kind of container that manages
// relationships between chain ID and CG Model.
template<typename realT>
class CGModelGenerator final
{
  public:
    typedef realT                           real_type;
    typedef CGChain<real_type>              cg_chain_type;
    typedef PDBChain<real_type>             pdb_chain_type;
    typedef CGModelGeneratorBase<real_type> generator_base;

  public:

    CGModelGenerator(): offset_(0){}
    CGModelGenerator(const std::size_t offset): offset_(offset){}
    CGModelGenerator(const CGModelGenerator&) = default;
    CGModelGenerator(CGModelGenerator&&)      = default;
    CGModelGenerator& operator=(const CGModelGenerator&) = default;
    CGModelGenerator& operator=(CGModelGenerator&&)      = default;

    template<typename ... Ts>
    void emplace_back(Ts&& ... vs)
    {
        this->generators_.emplace_back(std::forward<Ts>(vs)...);
    }

    std::vector<cg_chain_type>
    generate(const std::vector<pdb_chain_type>& pdbs) const
    {
        std::vector<cg_chain_type> cg_chains;

        std::size_t offset = this->offset_;
        for(const auto& id_gen : this->generators_)
        {
            const char  chain_id  = id_gen.first;
            const auto& generator = id_gen.second;

            const auto pdb_chain = std::find_if(pdbs.begin(), pdbs.end(),
                [=](const pdb_chain_type& pdbchn) noexcept -> bool {
                    return pdbchn.chain_id() == chain_id;
                });
            if(pdb_chain == pdbs.end())
            {
                mjolnir::throw_exception<std::runtime_error>("jarngreipr::"
                    "CGModelGenerator: expected chain ", chain_id, " is not "
                    "provided");
            }

            auto cg_chain = generator->generate(*pdb_chain, offset);
            offset += cg_chain.size();
            cg_chains.push_back(std::move(cg_chain));
        }
        return cg_chains;
    }

  private:

    std::size_t offset_;
    std::vector<
        // pairof chain-ID & ModelGenerator
        std::pair<char, std::shared_ptr<generator_base>>
        > generators_;
};

} // jarngreipr
#endif // JARNGREIPR_MODEL_CG_MODEL_GENERATOR_HPP
