#ifndef JARNGREIPR_CLEMENTI_GO
#define JARNGREIPR_CLEMENTI_GO
#include <jarngreipr/forcefield/ForceFieldGenerator.hpp>
#include <jarngreipr/geometry/distance.hpp>
#include <jarngreipr/geometry/angle.hpp>
#include <jarngreipr/geometry/dihedral.hpp>
#include <iterator>
#include <iostream>
#include <vector>

namespace jarngreipr
{

template<typename realT>
class ClementiGo final : public IntraChainForceFieldGenerator<realT>
{
  public:
    typedef IntraChainForceFieldGenerator<realT> base_type;
    typedef typename base_type::real_type real_type;
    typedef typename base_type::bead_type bead_type;

  public:

    ClementiGo()
        : contact_threshold_(6.5), k_bond_length_(100.0), k_bond_angle_(20.0),
          k_dihedral_angle_1_(1.0), k_dihedral_angle_3_(0.5), k_intra_go_(0.3)
    {}
    ~ClementiGo() override = default;

    explicit ClementiGo(const real_type th)
        : contact_threshold_(th), k_bond_length_(100.0), k_bond_angle_(20.0),
          k_dihedral_angle_1_(1.0), k_dihedral_angle_3_(0.5), k_intra_go_(0.3)
    {}

    ClementiGo(const real_type threshold,
        const real_type k_bl,   const real_type k_ba,
        const real_type k_dih1, const real_type k_dih3, const real_type k_igo)
        : contact_threshold_(threshold), k_bond_length_(k_bl), k_bond_angle_(k_ba),
          k_dihedral_angle_1_(k_dih1), k_dihedral_angle_3_(k_dih3),
          k_intra_go_(k_igo)
    {}

    // generate parameters and write out to `ostrm`.
    void generate(toml::Table& out,
        const std::vector<std::unique_ptr<bead_type>>& chain) const override;

    bool check_beads_kind(
        const std::vector<std::unique_ptr<bead_type>>& chain) const override;

    real_type& contact_threshold()       noexcept {return contact_threshold_;}
    real_type  contact_threshold() const noexcept {return contact_threshold_;}

  private:

    real_type min_distance_sq(
            const typename bead_type::container_type& lhs,
            const typename bead_type::container_type& rhs) const
    {
        real_type dist2 = std::numeric_limits<real_type>::max();
        for(const auto& l : lhs)
        {
            if(l.element == " H")
            {
                std::cerr << "ignore atom: " << l << std::endl;
                continue;
            }
            for(const auto& r : rhs)
            {
                if(r.element == " H")
                {
                    std::cerr << "ignore atom: " << r << std::endl;
                    continue;
                }
                const auto d2 = distance(l.position, r.position);
                dist2 = std::min(dist2, d2);
            }
        }
        return dist2;
    }

  private:
    real_type contact_threshold_;
    real_type k_bond_length_;
    real_type k_bond_angle_;
    real_type k_dihedral_angle_1_;
    real_type k_dihedral_angle_3_;
    real_type k_intra_go_;
};

template<typename realT>
void ClementiGo<realT>::generate(toml::Table& ff,
        const std::vector<std::unique_ptr<bead_type>>& chain) const
{
    if(!this->check_beads_kind(chain))
    {
        std::cerr << "ClementiGo: stop generating forcefield..." << std::endl;
        return ;
    }

    if(ff.count("local") == 0)
    {
        ff["local"] = toml::Array();
    }

    /* bond-length */ {
        toml::Table bond_length;
        bond_length["interaction"] = toml::String("BondLength");
        bond_length["potential"]   = toml::String("Harmonic");
        bond_length["topology"]    = toml::String("bond");

        toml::Array params;
        for(std::size_t i=0, sz = chain.size() - 1; i<sz; ++i)
        {
            const auto& bead1 = chain.at(i);
            const auto& bead2 = chain.at(i+1);
            const std::size_t i1 = bead1->index();
            const std::size_t i2 = bead2->index();
            toml::Table para;
            para["indices"] = toml::value{i1, i2};
            para["native"]  = distance(bead1->position(), bead2->position());
            para["k"]       = this->k_bond_length_;
            params.push_back(std::move(para));
        }
        bond_length["parameters"] = std::move(params);
        ff["local"].cast<toml::value_t::Array>().push_back(std::move(bond_length));
    }
    /* bond-angle */{
        toml::Table bond_angle;
        bond_angle["interaction"] = toml::String("BondAngle");
        bond_angle["potential"]   = toml::String("Harmonic");
        bond_angle["topology"]    = toml::String("none");

        toml::Array params;
        for(std::size_t i=0, sz = chain.size() - 2; i<sz; ++i)
        {
            const auto& bead1 = chain.at(i);
            const auto& bead2 = chain.at(i+1);
            const auto& bead3 = chain.at(i+2);
            const std::size_t i1 = bead1->index();
            const std::size_t i2 = bead2->index();
            const std::size_t i3 = bead3->index();

            toml::Table para;
            para["indices"] = toml::value{i1, i2, i3};
            para["native"]  = angle(bead1->position(), bead2->position(),
                                    bead3->position());
            para["k"]  = this->k_bond_angle_;
            params.push_back(std::move(para));
        }
        bond_angle["parameters"] = std::move(params);
        ff["local"].cast<toml::value_t::Array>().push_back(std::move(bond_angle));
    }
    /* dihedral-angle */{
        toml::Table dihd_angle;
        dihd_angle["interaction"] = toml::String("DihedralAngle");
        dihd_angle["potential"]   = toml::String("ClementiDihedral");
        dihd_angle["topology"]    = toml::String("none");

        toml::Array params;
        for(std::size_t i=0, sz = chain.size() - 3; i<sz; ++i)
        {
            const auto& bead1 = chain.at(i);
            const auto& bead2 = chain.at(i+1);
            const auto& bead3 = chain.at(i+2);
            const auto& bead4 = chain.at(i+3);
            const std::size_t i1 = bead1->index();
            const std::size_t i2 = bead2->index();
            const std::size_t i3 = bead3->index();
            const std::size_t i4 = bead4->index();

            toml::Table para;
            para["indices"] = toml::value{i1, i2, i3, i4};
            para["native"]  = dihedral_angle(bead1->position(),
                    bead2->position(), bead3->position(), bead4->position());
            para["k1"]  = this->k_dihedral_angle_1_;
            para["k3"]  = this->k_dihedral_angle_3_;
            params.push_back(std::move(para));
        }
        dihd_angle["parameters"] = std::move(params);
        ff["local"].cast<toml::value_t::Array>().push_back(std::move(dihd_angle));
    }

    const real_type th2 = this->contact_threshold_ * this->contact_threshold_;
    /* intra-chain-go-contacts */{
        toml::Table go_contact;
        go_contact["interaction"] = toml::String("BondLength");
        go_contact["potential"]   = toml::String("Go1012Contact");
        go_contact["topology"]    = toml::String("contact");

        toml::Array params;
        for(std::size_t i=0, sz_i = chain.size()-4; i<sz_i; ++i)
        {
            const auto& group1 = chain.at(i)->atoms();
            for(std::size_t j=i+4, sz_j = chain.size(); j<sz_j; ++j)
            {
                const auto& group2 = chain.at(j)->atoms();
                if(this->min_distance_sq(group1, group2))
                {
                    const auto& bead1 = chain.at(i);
                    const auto& bead2 = chain.at(j);
                    const std::size_t i1 = bead1->index();
                    const std::size_t i2 = bead2->index();

                    toml::Table para;
                    para["indices"] = toml::value{i1, i2};
                    para["native"]  = distance(bead1->position(), bead2->position());
                    para["k"]       = this->k_intra_go_;
                    params.push_back(std::move(para));
                }
            }
        }
        go_contact["parameters"]  = std::move(params);
        ff["local"].cast<toml::value_t::Array>().push_back(std::move(go_contact));
    }
    return;
}

template<typename realT>
bool ClementiGo<realT>::check_beads_kind(
        const std::vector<std::unique_ptr<bead_type>>& chain) const
{
    for(const auto& bead : chain)
    {
        if(bead->kind() != "CarbonAlpha")
        {
            std::cerr << "ClementiGo: invalid coarse-grained bead kind: "
                      << bead->kind() << '\n';
            std::cerr << "it allows only CarbonAlpha beads.\n";
            return false;
        }
    }
    return true;
}

}//jarngreipr
#endif /* JARNGREIPR_CLEMENTI_GO */
