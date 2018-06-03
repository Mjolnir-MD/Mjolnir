#ifndef JARNGREIPR_MODEL_GRAINER_HPP
#define JARNGREIPR_MODEL_GRAINER_HPP
#include <jarngreipr/model/Bead.hpp>
#include <jarngreipr/model/CGChain.hpp>
#include <jarngreipr/pdb/PDBChain.hpp>

namespace jarngreipr
{

template<typename realT>
class GrainerBase
{
  public:
    typedef realT           real_type;
    typedef Bead<real_type> bead_type;

  public:

    virtual ~GrainerBase() = default;

    virtual CGChain<realT>
    grain(const PDBChain<realT>& pdb, const std::size_t offset) const = 0;
};

template<typename realT>
class Grainer
{
  public:
    typedef realT           real_type;
    typedef Bead<real_type> bead_type;

  public:

    Grainer() = default;

    void push_back(std::shared_ptr<GrainerBase<real_type>>&& item)
    {
        this->grainers_.push_back(std::move(item));
        return;
    }

    std::vector<CGChain<real_type>>
    grain(const std::vector<PDBChain<realT>>& pdbs) const
    {
        if(pdbs.size() != grainers_.size())
        {
            mjolnir::throw_exception<std::runtime_error>("size of pdb chains(",
                pdbs.size(), ") differ from the number of model specifications",
                this->grainers_.size());
        }

        std::vector<CGChain<realT>> cg_chains;
        std::size_t offset = 0;
        for(std::size_t i=0; i<pdbs.size(); ++i)
        {
            auto cg_chain = this->grainers_.at(i)->grain(pdbs.at(i), offset);
            offset += cg_chain.size();
            cg_chains.push_back(std::move(cg_chain));
        }

        return cg_chains;
    }

  private:

    std::vector<std::shared_ptr<GrainerBase<real_type>>> grainers_;
};


} // jarngreipr
#endif // JARNGREIPR_MODEL_GRAINER
