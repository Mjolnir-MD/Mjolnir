#ifndef JARNGREIPR_MODEL_CARBON_ALPHA_HPP
#define JARNGREIPR_MODEL_CARBON_ALPHA_HPP
#include <jarngreipr/model/Bead.hpp>
#include <jarngreipr/model/CGModelGenerator.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <algorithm>
#include <stdexcept>
#include <memory>
#include <string>

namespace jarngreipr
{

/*! @brief carbon alpha 1 beads per amino acid model */
template<typename realT>
class CarbonAlpha final : public Bead<realT>
{
  public:
    typedef Bead<realT> base_type;
    typedef typename base_type::real_type       real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::atom_type       atom_type;
    typedef typename base_type::container_type  container_type;

  public:

    CarbonAlpha(std::size_t idx, container_type atoms, std::string name)
        : base_type(idx, std::move(atoms), std::move(name))
    {
        if(!this->atoms_.empty())
        {
            const auto is_ca =
                [](const atom_type& a){return a.atom_name==" CA ";};
            const std::size_t num_ca = std::count_if(
                this->atoms_.cbegin(), this->atoms_.cend(), is_ca);
            if(num_ca == 0)
            {
                mjolnir::throw_exception<std::runtime_error>("jarngreipr::"
                    "model::CarbonAlpha: no c-alpha atom in this residue: \n",
                    this->atoms_.front());
            }
            if(num_ca > 1)
            {
                mjolnir::throw_exception<std::runtime_error>("jarngreipr::"
                    "model::CarbonAlpha: multiple c-alpha in this residue: \n",
                    this->atoms_.front());
            }
            this->position_ = std::find_if(
                this->atoms_.cbegin(), this->atoms_.cend(), is_ca)->position;
        }
    }
    ~CarbonAlpha() override = default;

    CarbonAlpha(const CarbonAlpha&) = default;
    CarbonAlpha(CarbonAlpha&&)      = default;
    CarbonAlpha& operator=(const CarbonAlpha&) = default;
    CarbonAlpha& operator=(CarbonAlpha&&)      = default;

    std::string attribute(const std::string& n) const override {return "";}
    std::string kind() const override {return "CarbonAlpha";}

    coordinate_type position() const override {return this->position_;}

  private:

    coordinate_type position_;
};

template<typename realT>
class CarbonAlphaGenerator final : public CGModelGeneratorBase<realT>
{
  public:
    typedef CGModelGeneratorBase<realT>        base_type;
    typedef typename base_type::real_type      real_type;
    typedef typename base_type::cg_chain_type  cg_chain_type;
    typedef typename base_type::pdb_chain_type pdb_chain_type;

  public:

    CarbonAlphaGenerator() = default;
    ~CarbonAlphaGenerator() override = default;

    cg_chain_type
    generate(const pdb_chain_type& pdb, const std::size_t offset) const override
    {
        CGChain<realT> retval(std::string(1, pdb.chain_id()));
        for(std::size_t i=0; i<pdb.residues_size(); ++i)
        {
            const auto res = pdb.residue_at(i);
            std::vector<PDBAtom<realT>> atoms(res.begin(), res.end());
            const auto name = atoms.front().residue_name;
            retval.push_back(std::make_shared<CarbonAlpha<realT>>(
                        i + offset, std::move(atoms), name));
        }
        return retval;
    }
};

} // jarngreipr
#endif /*JARNGREIPR_CARBON_ALPHA*/
