#ifndef JARNGREIPR_FORCEFIELD_EXCLUDED_VOLUME
#define JARNGREIPR_FORCEFIELD_EXCLUDED_VOLUME
#include <jarngreipr/forcefield/ForceFieldGenerator.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <extlib/toml/toml.hpp>

namespace jarngreipr
{

template<typename realT>
class ExcludedVolume final : public ForceFieldGenerator<realT>
{
  public:
    typedef InterChainForceFieldGenerator<realT> base_type;
    typedef typename base_type::real_type real_type;
    typedef typename base_type::bead_type bead_type;
    typedef typename base_type::bead_ptr  bead_ptr;

  public:

    ExcludedVolume(const toml::Table& para) : parameters_(para){}
    ~ExcludedVolume() override = default;

    void generate(toml::Table& out,
        const std::vector<std::vector<bead_ptr>>& chains) const;

    bool check_beads_kind(
        const std::vector<std::vector<bead_type>>& chain) const {return true;}

  private:

    toml::Table parameters_;
};

template<typename realT>
void ExcludedVolume<realT>::generate(toml::Table& ff,
        const std::vector<std::vector<bead_ptr>>& chains) const
{
    if(ff.count("global") == 0)
    {
        ff["global"] = toml::Array();
    }

    toml::Table exv;
    exv["interaction"]      = toml::String("Distance");
    exv["potential"]        = toml::String("ExcludedVolume");
    exv["ignored_chain"]    = toml::String("Nothing");
    exv["ignored_bonds"]    = 3;
    exv["ignored_contacts"] = 1;

    toml::Table partition;
    partition["type"]   = toml::String("CellList");
    partition["margin"] = 1.0;
    exv["spatial_partition"] = partition;

    toml::Array params;
    for(const auto& chain : chains)
    {
        for(const auto& bead : chain)
        {
            toml::Table para;
            para["index"] = bead->index();
            para["sigma"] = toml::get<toml::Float>(mjolnir::toml_value_at(
                this->parameters_, bead->name(), "jarngreipr::ExcludedVolume"));
            params.push_back(para);
        }
    }
    exv["epsilon"]    = 0.2;
    exv["parameters"] = toml::value(std::move(params));

    ff["global"].cast<toml::value_t::Array>().push_back(exv);
    return;
}

} // jarngreipr
#endif// JARNGREIPR_FORCEFIELD_EXCLUDED_VOLUME
