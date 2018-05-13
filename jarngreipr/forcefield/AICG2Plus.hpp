#ifndef JARNGREIPR_FORCEFIELD_AICG2_PLUS_H
#define JARNGREIPR_FORCEFIELD_AICG2_PLUS_H
#include <jarngreipr/forcefield/ForceFieldGenerator.hpp>
#include <jarngreipr/geometry/distance.hpp>
#include <jarngreipr/geometry/angle.hpp>
#include <jarngreipr/geometry/dihedral.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <iterator>
#include <algorithm>
#include <iostream>
#include <vector>
#include <limits>

namespace jarngreipr
{

template<typename realT>
class AICG2Plus final : public IntraChainForceFieldGenerator<realT>
{
  public:
    typedef IntraChainForceFieldGenerator<realT> base_type;
    typedef typename base_type::real_type real_type;
    typedef typename base_type::bead_type bead_type;

  public:

    AICG2Plus(const toml::Table& para, const std::vector<std::size_t>& flex);
    ~AICG2Plus() override = default;

    // generate parameters and write out to `ostrm`.
    void generate(toml::Table& out,
        const std::vector<std::unique_ptr<bead_type>>& chain) const override;

    bool check_beads_kind(
        const std::vector<std::unique_ptr<bead_type>>& chain) const override;

  private:

    bool is_flexible_region(const std::size_t bead_idx) const
    {
        return std::binary_search(this->flexible_beads_.begin(),
                                  this->flexible_beads_.end(), bead_idx);
    }

    bool is_backbone(const PDBAtom<realT>& atom) const
    {
        return (atom.atom_name == " N  " || atom.atom_name == " C  " ||
                atom.atom_name == " O  " || atom.atom_name == " OXT" ||
                atom.atom_name == " CA ");
    }
    bool is_sidechain(const PDBAtom<realT>& atom) const
    {
        return !is_backbone(atom) &&
            atom.atom_name.at(0) != 'H' && atom.atom_name.at(1) != 'H';
    }

    bool is_donor(const PDBAtom<realT>& atom) const
    {
        return (atom.atom_name.at(1) == 'N' ||
                (atom.residue_name == "SER" && atom.atom_name == " OG ") ||
                (atom.residue_name == "THR" && atom.atom_name == " OG1") ||
                (atom.residue_name == "TYR" && atom.atom_name == " OH ") ||
                (atom.residue_name == "CYS" && atom.atom_name.at(1) == 'S'));
    }
    bool is_acceptor(const PDBAtom<realT>& atom) const
    {
        return (atom.atom_name.at(1) == 'O' || atom.atom_name.at(1) == 'S');
    }
    bool is_cation(const PDBAtom<realT>& atom) const
    {
	    return ((atom.residue_name == "ARG" && atom.atom_name == " NH1") ||
                (atom.residue_name == "ARG" && atom.atom_name == " NH2") ||
                (atom.residue_name == "LYS" && atom.atom_name == " NZ "));
    }
    bool is_anion(const PDBAtom<realT>& atom) const
    {
		return ((atom.residue_name == "GLU" && atom.atom_name == " OE1") ||
                (atom.residue_name == "GLU" && atom.atom_name == " OE2") ||
                (atom.residue_name == "ASP" && atom.atom_name == " OD1") ||
                (atom.residue_name == "ASP" && atom.atom_name == " OD2"));
    }
    bool is_carbon(const PDBAtom<realT>& atom) const
    {
        return (atom.atom_name.at(1) == 'C');
    }

    real_type calc_contact_coef(const std::unique_ptr<bead_type>& bead1,
                                const std::unique_ptr<bead_type>& bead2) const;

    real_type min_distance_sq(const std::unique_ptr<bead_type>& bead1,
                              const std::unique_ptr<bead_type>& bead2) const
    {
        real_type min_dist = std::numeric_limits<real_type>::max();
        for(const auto& atom1 : bead1->atoms())
        {
            if(atom1.element == " H") {continue;}
            for(const auto& atom2 : bead2->atoms())
            {
                if(atom2.element == " H") {continue;}
                const real_type dist = distance_sq(atom1.position, atom2.position);
                if(dist < min_dist) {min_dist = dist;}
            }
        }
        return min_dist;
    }

    real_type limit_energy(const real_type val) const
    {
        if(val >= this->e_max_) {return e_max_;}
        if(val <= this->e_min_) {return e_min_;}
        return val;
    }

  private:

    // cache
    toml::Float e_min_;
    toml::Float e_max_;

    toml::Float go_contact_threshold_;
    toml::Float atom_contact_cutoff_;
    toml::Float hydrogen_bond_cutoff_;
    toml::Float solt_bridge_cutoff_;

    toml::Float coef_13_;
    toml::Float coef_14_;
    toml::Float coef_go_;

    toml::Float bb_hydrogen_bond_;
    toml::Float bb_donor_acceptor_;
    toml::Float bb_carbon_contact_;
    toml::Float bb_other_contact_;
    toml::Float ss_hydrogen_bond_;
    toml::Float ss_donor_acceptor_;
    toml::Float ss_salty_bridge_;
    toml::Float ss_carbon_contact_;
    toml::Float ss_charge_contact_;
    toml::Float ss_other_contact_;
    toml::Float bs_hydrogen_bond_;
    toml::Float bs_donor_acceptor_;
    toml::Float bs_carbon_contact_;
    toml::Float bs_charge_contact_;
    toml::Float bs_other_contact_;
    toml::Float long_range_contact_;
    toml::Float offset_;

    toml::Table parameters_; // [AICG2Plus] section
    std::vector<std::size_t> flexible_beads_;
};

template<typename realT>
void AICG2Plus<realT>::generate(toml::Table& ff,
        const std::vector<std::unique_ptr<bead_type>>& chain) const
{
    if(!this->check_beads_kind(chain))
    {
        std::cerr << "AICG2+: Invalid Bead Kind. stop generating forcefield..."
                  << std::endl;
        return ;
    }

    if(ff.count("local") == 0)
    {
        ff["local"] = toml::Array();
    }

    /* bond-length */ {
        toml::Table bond_length;
        bond_length["interaction"] = toml::String("BondLength");
        bond_length["potential"  ] = toml::String("Harmonic");
        bond_length["topology"   ] = toml::String("bond");

        const auto& kbd = mjolnir::toml_value_at(
                this->parameters_, "cbd_aicg2", "jarngreipr::AICG2Plus");

        toml::Array params;
        for(std::size_t i=0, sz = chain.size() - 1; i<sz; ++i)
        {
            const auto& bead1 = chain.at(i);
            const auto& bead2 = chain.at(i+1);
            const std::size_t i1 = bead1->index();
            const std::size_t i2 = bead2->index();
            toml::Table para;
            para["indices"] = toml::value{i1, i2};
            para["eq"     ] = distance(bead1->position(), bead2->position());
            para["k"      ] = kbd;
            params.push_back(std::move(para));
        }
        bond_length["parameters"] = std::move(params);
        ff["local"].cast<toml::value_t::Array>().push_back(std::move(bond_length));
    }
    /* bond-angle */{
        toml::Table bond_angle;
        bond_angle["interaction"] = toml::String("BondLength");
        bond_angle["potential"]   = toml::String("AICG2PlusAngle");
        bond_angle["topology"]    = toml::String("none");

        const auto& width = mjolnir::toml_value_at(
                this->parameters_, "wid_aicg13", "jarngreipr::AICG2Plus");

        toml::Array params;
        for(std::size_t i=0, sz = chain.size() - 2; i<sz; ++i)
        {
            const auto& bead1 = chain.at(i);
            const auto& bead2 = chain.at(i+1);
            const auto& bead3 = chain.at(i+2);
            const std::size_t i1 = bead1->index();
            const std::size_t i2 = bead2->index();
            const std::size_t i3 = bead3->index();

            // if the beads contains flexible region, remove 1-3 contact.
            if(is_flexible_region(i1) || is_flexible_region(i2) ||
               is_flexible_region(i3)) {continue;}

            toml::Table para;
            para["indices"] = toml::value{i1, i3};
            para["eq"     ] = distance(bead1->position(), bead3->position());
            para["w"      ] = width;
            para["epsilon"] = this->coef_13_ * calc_contact_coef(bead1, bead3);
            params.push_back(std::move(para));
        }
        bond_angle["parameters"] = std::move(params);
        ff["local"].cast<toml::value_t::Array>().push_back(std::move(bond_angle));
    }
    /* flexible-local-angle */{
        toml::Table flp_angle;
        flp_angle["interaction"] = toml::String("BondAngle");
        flp_angle["potential"  ] = toml::String("FlexibleLocalAngle");
        flp_angle["topology"   ] = toml::String("none");

        const auto& flp = mjolnir::toml_value_at(
            this->parameters_, "flexible_local", "jarngreipr::AICG2Plus"
            ).template cast<toml::value_t::Table>();

        const auto& term1 = mjolnir::toml_value_at(
            flp, "angle_term1", "jarngreipr::AICG2Plus"
            ).template cast<toml::value_t::Table>();
        const auto& term2 = mjolnir::toml_value_at(
            flp, "angle_term2", "jarngreipr::AICG2Plus"
            ).template cast<toml::value_t::Table>();

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
            para["k"      ] = mjolnir::toml_value_at(flp, "k_angle",
                                                     "jarngreipr::AICG2Plus");
            para["term1"  ] = mjolnir::toml_value_at(term1, bead2->name(),
                                                     "jarngreipr::AICG2Plus");
            para["term2"  ] = mjolnir::toml_value_at(term2, bead2->name(),
                                                     "jarngreipr::AICG2Plus");
            params.push_back(std::move(para));
        }
        flp_angle["parameters"] = std::move(params);
        ff["local"].cast<toml::value_t::Array>().push_back(std::move(flp_angle));
    }
    /* dihedral-angle */{
        toml::Table dihd_angle;
        dihd_angle["interaction"] = toml::String("DihedralAngle");
        dihd_angle["potential"  ] = toml::String("AICG2PlusDihedral");
        dihd_angle["topology"   ] = toml::String("none");

        const auto& width = mjolnir::toml_value_at(
            this->parameters_, "wid_dih", "jarngreipr::AICG2Plus");

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

            // if the beads contains flexible region, remove 1-4 contact.
            if(is_flexible_region(i1) || is_flexible_region(i2) ||
               is_flexible_region(i3) || is_flexible_region(i4)) {continue;}

            toml::Table para;
            para["indices"] = toml::value{i1, i2, i3, i4};
            para["eq"     ] = dihedral_angle(bead1->position(),
                    bead2->position(), bead3->position(), bead4->position());
            para["w"      ] = width;
            para["epsilon"] = this->coef_14_ * calc_contact_coef(bead1, bead4);
            params.push_back(std::move(para));
        }
        dihd_angle["parameters"] = std::move(params);
        ff["local"].cast<toml::value_t::Array>().push_back(std::move(dihd_angle));
    }
    /* flexible-dihedral-angle */{
        toml::Table dihd_angle;
        dihd_angle["interaction"] = toml::String("DihedralAngle");
        dihd_angle["potential"  ] = toml::String("FlexibleLocalDihedral");
        dihd_angle["topology"   ] = toml::String("none");

        const auto& flp = mjolnir::toml_value_at(
            this->parameters_, "flexible_local", "jarngreipr::AICG2Plus"
            ).template cast<toml::value_t::Table>();

        const auto& term = mjolnir::toml_value_at(
            flp, "dihedral_term", "jarngreipr::AICG2Plus"
            ).template cast<toml::value_t::Table>();

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

            // like "ALA-PHE" or something like that
            const std::string ident =
                bead2->name() + std::string("-") + bead3->name();

            toml::Table para;
            para["indices"] = toml::value{i1, i2, i3, i4};
            para["k"      ] = mjolnir::toml_value_at(flp, "k_dihedral",
                                                     "jarngreipr::AICG2Plus");
            para["term"   ] = mjolnir::toml_value_at(term, ident,
                                                     "jarngreipr::AICG2Plus");
            params.push_back(std::move(para));
        }
        dihd_angle["parameters"] = std::move(params);
        ff["local"].cast<toml::value_t::Array>().push_back(std::move(dihd_angle));
    }

    const real_type th2 = this->go_contact_threshold_ * this->go_contact_threshold_;
    /* intra-chain-go-contacts */{
        toml::Table go_contact;
        go_contact["interaction"] = toml::String("BondLength");
        go_contact["potential"  ] = toml::String("Go1012Contact");
        go_contact["topology"   ] = toml::String("contact");

        toml::Array params;
        for(std::size_t i=0, sz_i = chain.size()-4; i<sz_i; ++i)
        {
            for(std::size_t j=i+4, sz_j = chain.size(); j<sz_j; ++j)
            {
                if(this->min_distance_sq(chain.at(i), chain.at(j)) < th2)
                {
                    const auto& bead1 = chain.at(i);
                    const auto& bead2 = chain.at(j);
                    const std::size_t i1 = bead1->index();
                    const std::size_t i2 = bead2->index();

                    // if one of the bead is flexible region, continue.
                    if(is_flexible_region(i1) || is_flexible_region(i2))
                    {continue;}

                    toml::Table para;
                    para["indices"] = toml::value{i1, i2};
                    para["eq"     ] = distance(bead1->position(), bead2->position());
                    para["k"      ] = -this->coef_go_ * calc_contact_coef(bead1, bead2);
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
typename AICG2Plus<realT>::real_type
AICG2Plus<realT>::calc_contact_coef(
        const std::unique_ptr<bead_type>& bead1,
        const std::unique_ptr<bead_type>& bead2) const
{
    std::size_t num_bb_hb = 0; // hydrogen bond
    std::size_t num_bb_da = 0; // donor-acceptor
    std::size_t num_bb_cx = 0; // carbon-X contact
    std::size_t num_bb_oc = 0; // other contact

    std::size_t num_ss_hb = 0; // hydrogen bond
    std::size_t num_ss_sb = 0; // salty bridge
    std::size_t num_ss_da = 0; // donor-acceptor
    std::size_t num_ss_cc = 0; // charge contact
    std::size_t num_ss_cx = 0; // carbon-X contact
    std::size_t num_ss_oc = 0; // other contact

    std::size_t num_bs_hb = 0; // hydrogen bond
    std::size_t num_bs_da = 0; // donor-acceptor
    std::size_t num_bs_cc = 0; // charge contact
    std::size_t num_bs_cx = 0; // carbon-X contact
    std::size_t num_bs_oc = 0; // other contact

    std::size_t num_long  = 0; // long range contact

    for(const auto& atom1 : bead1->atoms())
    {
    for(const auto& atom2 : bead2->atoms())
    {
        const auto dist = distance(atom1.position, atom2.position);
        if(dist < go_contact_threshold_) {num_long += 1;}
        if(dist >= atom_contact_cutoff_) {continue;}

        num_long  -= 1; // it is short range contact. reduce the contribution.
        if(is_backbone(atom1) && is_backbone(atom2))
        {
            if((is_acceptor(atom1) && is_donor(atom2)) ||
               (is_acceptor(atom2) && is_donor(atom1)))
            {
                if(dist < hydrogen_bond_cutoff_)
                {
                    num_bb_hb += 1;
                }
                else // no hydrogen bond, just donor-acceptor
                {
                    num_bb_da += 1;
                }
            }
            else if(is_carbon(atom1) || is_carbon(atom2))
            {
                num_bb_cx += 1;
            }
            else // other atomic contacts
            {
                num_bb_oc += 1;
            }
        }
        else if(is_sidechain(atom1) && is_sidechain(atom2))
        {
            if((is_acceptor(atom1) && is_donor(atom2)) ||
               (is_acceptor(atom2) && is_donor(atom1)))
            {
                if((is_cation(atom1) && is_anion(atom2)) ||
                   (is_cation(atom2) && is_anion(atom1)))
                {
                    if(dist < solt_bridge_cutoff_)
                    {
                        num_ss_sb += 1;
                    }
                    else // not a solt bridge, charge contact.
                    {
                        num_ss_cc += 1;
                    }
                }
                else if(dist < hydrogen_bond_cutoff_)
                {
                    num_ss_hb += 1;
                }
                else if(is_cation(atom1) || is_anion(atom2) ||
                        is_cation(atom2) || is_anion(atom1))
                {
                    // one of the two has a charge, charge contact.
                    num_ss_cc += 1;
                }
                else // just a donor-acceptor.
                {
                    num_ss_da += 1;
                }
            }
            else if(is_cation(atom1) || is_anion(atom2) ||
                    is_cation(atom2) || is_anion(atom1))
            {
                num_ss_cc += 1; // charge contact.
            }
            else if(is_carbon(atom1) || is_carbon(atom2))
            {
                num_ss_cx += 1; // carbon-X
            }
            else // other atomic contacts
            {
                num_ss_oc += 1;
            }
        }
        else if((is_sidechain(atom1) && is_backbone(atom2)) ||
                (is_sidechain(atom2) && is_backbone(atom1)))
        {
            if((is_acceptor(atom1) && is_donor(atom2)) ||
               (is_acceptor(atom2) && is_donor(atom1)))
            {
                if(dist < hydrogen_bond_cutoff_)
                {
                    num_bs_hb += 1; // hydrogen bond formed.
                }
                else if(is_cation(atom1) || is_anion(atom2) ||
                        is_cation(atom2) || is_anion(atom1))
                {
                    num_bs_cc += 1; // one of two has a charge.
                }
                else
                {
                    num_bs_da += 1; // just donor-acceptor.
                }
            }
            else if(is_cation(atom1) || is_anion(atom2) ||
                    is_cation(atom2) || is_anion(atom1))
            {
                num_bs_cc += 1; // one of two has a charge.
            }
            else if(is_carbon(atom1) || is_carbon(atom2))
            {
                num_bs_cx += 1; // one of two is a carbon
            }
            else // other kind of atomic contacts...
            {
                num_bs_oc += 1;
            }
        }
    } // atom2
    } // atom1

    // XXX pair of residues cannot form more than one salty bridge.
    // set the number of salty bridge to 1, and the rests are counted as
    // sidechain-sidechain charge contact(backbone never forms salty bridge).
    if(num_ss_sb > 1)
    {
        num_ss_cc += (num_ss_sb - 1);
        num_ss_sb = 1;
    }

    // calculate weighted sum.
    real_type e_tot = offset_;
    e_tot += bb_hydrogen_bond_   * num_bb_hb;
    e_tot += bb_donor_acceptor_  * num_bb_da;
    e_tot += bb_carbon_contact_  * num_bb_cx;
    e_tot += bb_other_contact_   * num_bb_oc;
    e_tot += ss_hydrogen_bond_   * num_ss_hb;
    e_tot += ss_donor_acceptor_  * num_ss_da;
    e_tot += ss_salty_bridge_    * num_ss_sb;
    e_tot += ss_carbon_contact_  * num_ss_cx;
    e_tot += ss_charge_contact_  * num_ss_cc;
    e_tot += ss_other_contact_   * num_ss_oc;
    e_tot += bs_hydrogen_bond_   * num_bs_hb;
    e_tot += bs_donor_acceptor_  * num_bs_da;
    e_tot += bs_carbon_contact_  * num_bs_cx;
    e_tot += bs_charge_contact_  * num_bs_cc;
    e_tot += bs_other_contact_   * num_bs_oc;
    e_tot += long_range_contact_ * num_long;

    return this->limit_energy(e_tot);
}

template<typename realT>
bool AICG2Plus<realT>::check_beads_kind(
        const std::vector<std::unique_ptr<bead_type>>& chain) const
{
    // TODO? its good to check all the beads has appropreate atoms
    for(const auto& bead : chain)
    {
        if(bead->kind() != "CarbonAlpha")
        {
            std::cerr << "AICG2Plus: invalid coarse-grained bead kind: "
                      << bead->kind() << '\n';
            std::cerr << "it allows only CarbonAlpha beads.\n";
            return false;
        }
    }
    return true;
}

template<typename realT>
AICG2Plus<realT>::AICG2Plus(
        const toml::Table& para, const std::vector<std::size_t>& flex)
    : parameters_(para), flexible_beads_(flex)
{
    // to use std::binary_search later
    std::sort(this->flexible_beads_.begin(), this->flexible_beads_.end());

    // read and cache the parameters to be used...

    this->e_min_ = toml::get<toml::Float>(mjolnir::toml_value_at(
            this->parameters_, "ecut_up_aicg2", "[AICG2+]"));
    this->e_max_ = toml::get<toml::Float>(mjolnir::toml_value_at(
            this->parameters_, "ecut_low_aicg2",  "[AICG2+]"));

    this->go_contact_threshold_ = toml::get<toml::Float>(
        mjolnir::toml_value_at(this->parameters_,
            "go_contact_threshold", "[AICG2+]"));
    this->atom_contact_cutoff_ = toml::get<toml::Float>(
        mjolnir::toml_value_at(this->parameters_,
            "atom_contact_cutoff", "[AICG2+]"));
    this->solt_bridge_cutoff_ = toml::get<toml::Float>(
        mjolnir::toml_value_at(this->parameters_,
            "solt_bridge_cutoff", "[AICG2+]"));
    this->hydrogen_bond_cutoff_ = toml::get<toml::Float>(
        mjolnir::toml_value_at(this->parameters_,
            "hydrogen_bond_cutoff", "[AICG2+]"));

    this->coef_13_ = toml::get<toml::Float>(mjolnir::toml_value_at(
                this->parameters_, "caicg2plus_13", "[AICG2+]"));
    this->coef_14_ = toml::get<toml::Float>(mjolnir::toml_value_at(
                this->parameters_, "caicg2plus_14", "[AICG2+]"));
    this->coef_go_ = toml::get<toml::Float>(mjolnir::toml_value_at(
                this->parameters_, "caicg2plus_nloc", "[AICG2+]"));

    const auto& contact_energy_coef = mjolnir::toml_value_at(this->parameters_,
            "contact_energy_coefficients", "jarngreipr::AICG2Plus"
            ).template cast<toml::value_t::Table>();

    this->bb_hydrogen_bond_  = toml::get<toml::Float>(
        mjolnir::toml_value_at(contact_energy_coef, "backbone_hydrogen_bond",
            "[AICG2+.contact_energy_coefficients]"));
    this->bb_donor_acceptor_ = toml::get<toml::Float>(
        mjolnir::toml_value_at(contact_energy_coef, "backbone_donor_acceptor",
            "[AICG2+.contact_energy_coefficients]"));
    this->bb_carbon_contact_ = toml::get<toml::Float>(
        mjolnir::toml_value_at(contact_energy_coef, "backbone_carbon_contact",
            "[AICG2+.contact_energy_coefficients]"));
    this->bb_other_contact_  = toml::get<toml::Float>(
        mjolnir::toml_value_at(contact_energy_coef, "backbone_contact",
            "[AICG2+.contact_energy_coefficients]"));

    this->ss_hydrogen_bond_  = toml::get<toml::Float>(
        mjolnir::toml_value_at(contact_energy_coef, "sidechain_hydrogen_bond",
            "[AICG2+.contact_energy_coefficients]"));
    this->ss_donor_acceptor_ = toml::get<toml::Float>(
        mjolnir::toml_value_at(contact_energy_coef, "sidechain_donor_acceptor",
            "[AICG2+.contact_energy_coefficients]"));
    this->ss_salty_bridge_   = toml::get<toml::Float>(
        mjolnir::toml_value_at(contact_energy_coef, "sidechain_salty_bridge",
            "[AICG2+.contact_energy_coefficients]"));
    this->ss_carbon_contact_ = toml::get<toml::Float>(
        mjolnir::toml_value_at(contact_energy_coef, "sidechain_carbon_contact",
            "[AICG2+.contact_energy_coefficients]"));
    this->ss_charge_contact_ = toml::get<toml::Float>(
        mjolnir::toml_value_at(contact_energy_coef, "sidechain_charge_contact",
            "[AICG2+.contact_energy_coefficients]"));
    this->ss_other_contact_  = toml::get<toml::Float>(
        mjolnir::toml_value_at(contact_energy_coef, "sidechain_contact",
            "[AICG2+.contact_energy_coefficients]"));

    this->bs_hydrogen_bond_  = toml::get<toml::Float>(
        mjolnir::toml_value_at(contact_energy_coef, "heterogeneous_hydrogen_bond",
            "[AICG2+.contacgy_coefficients]"));
    this->bs_donor_acceptor_ = toml::get<toml::Float>(
        mjolnir::toml_value_at(contact_energy_coef, "heterogeneous_donor_acceptor",
            "[AICG2+.contact_energy_coefficients]"));
    this->bs_carbon_contact_ = toml::get<toml::Float>(
        mjolnir::toml_value_at(contact_energy_coef, "heterogeneous_carbon_contact",
            "[AICG2+.contact_energy_coefficients]"));
    this->bs_charge_contact_ = toml::get<toml::Float>(
        mjolnir::toml_value_at(contact_energy_coef, "heterogeneous_charge_contact",
            "[AICG2+.contact_energy_coefficients]"));
    this->bs_other_contact_  = toml::get<toml::Float>(
        mjolnir::toml_value_at(contact_energy_coef, "heterogeneous_contact",
            "[AICG2+.contact_energy_coefficients]"));

    this->long_range_contact_ = toml::get<toml::Float>(
        mjolnir::toml_value_at(contact_energy_coef, "long_range_contact",
            "[AICG2+.contact_energy_coefficients]"));

    this->offset_ = toml::get<toml::Float>(mjolnir::toml_value_at(
        contact_energy_coef, "offset", "[AICG2+.contact_energy_coefficients]"));
}

}//jarngreipr
#endif /* JARNGREIPR_FORCEFIELD_AICG2_PLUS */
