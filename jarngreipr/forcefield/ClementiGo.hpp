#ifndef JARNGREIPR_CLEMENTI_GO
#define JARNGREIPR_CLEMENTI_GO
#include <jarngreipr/forcefield/ForceFieldGenerator.hpp>
#include <jarngreipr/geometry/distance.hpp>
#include <jarngreipr/geometry/angle.hpp>
#include <jarngreipr/geometry/dihedral.hpp>
#include <iterator>
#include <iostream>
#include <vector>

namespace mjolnir
{

template<typename coordT>
class ClementiGo final : public IntraChainForceFieldGenerator<coordT>
{
  public:
    typedef coordT coordinate_type;
    typedef IntraChainForceFieldGenerator<coordT> base_type;
    typedef typename base_type::atom_type         atom_type;
    typedef typename base_type::residue_type      residue_type;
    typedef typename base_type::chain_type        chain_type;
    typedef typename base_type::cg_chain_type     cg_chain_type;
    typedef typename base_type::connection_info   connection_info;
    typedef typename scalar_type_of<coordinate_type>::type real_type;

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

    ClementiGo(const real_type conth,
            const real_type k_bl, const real_type k_ba, const real_type k_dih1,
            const real_type k_dih3, const real_type k_igo)
        : contact_threshold_(conth), k_bond_length_(k_bl), k_bond_angle_(k_ba),
          k_dihedral_angle_1_(k_dih1), k_dihedral_angle_3_(k_dih3),
          k_intra_go_(k_igo)
    {}

    // generate parameters and write out to `ostrm`.
    connection_info
    generate(std::ostream& ostrm, const cg_chain_type& chain,
             const std::size_t offset) const override;

    bool check_beads_kind(const cg_chain_type& chain) const override;

    real_type& contact_threshold()       noexcept {return contact_threshold_;}
    real_type  contact_threshold() const noexcept {return contact_threshold_;}

  private:
    real_type contact_threshold_;
    real_type k_bond_length_;
    real_type k_bond_angle_;
    real_type k_dihedral_angle_1_;
    real_type k_dihedral_angle_3_;
    real_type k_intra_go_;
};

template<typename coordT>
typename ClementiGo<coordT>::connection_info
ClementiGo<coordT>::generate(std::ostream& ostrm,
        const cg_chain_type& chain, const std::size_t offset) const
{
    if(false == this->check_beads_kind(chain))
    {
        throw std::invalid_argument("jarngreipr::ClementiGo::generate: "
                "invalid bead kind appear in argument `chain`.");
    }
    connection_info connections;
    for(std::size_t i=0; i < chain.size(); ++i)
    {
        const std::size_t index = i + offset;
        connections[index].insert(index);
    }

    const real_type th2 = this->contact_threshold_ * this->contact_threshold_;

    /* bond-length */ {
        ostrm << "[[forcefields.local]]\n";
        ostrm << "interaction = \"BondLength\"\n";
        ostrm << "potential   = \"Harmonic\"\n";
        ostrm << "parameters  = [\n";
        for(std::size_t i=0, sz = chain.size() - 1; i<sz; ++i)
        {
            const std::size_t index1 = i + offset;
            const std::size_t index2 = i + offset + 1;

            connections[index1].insert(index2);
            connections[index2].insert(index1);

            ostrm << "{indices = [" << index1 << ", " << index2 << "], ";
            ostrm << "native = " << std::fixed << std::showpoint
                  << distance(chain.at(i), chain.at(i+1)) << ", ";
            ostrm << "k = " << std::fixed << std::showpoint
                  << this->k_bond_length_;
            ostrm << "},\n";
        }
        ostrm << "]\n";
    }
    /* bond-angle */{
        ostrm << "[[forcefields.local]]\n";
        ostrm << "interaction = \"BondAngle\"\n";
        ostrm << "potential   = \"Harmonic\"\n";
        ostrm << "parameters  = [\n";
        for(std::size_t i=0, sz = chain.size() - 2; i<sz; ++i)
        {
            const std::size_t index1 = i + offset;
            const std::size_t index2 = i + offset + 1;
            const std::size_t index3 = i + offset + 2;

            ostrm << "{indices = [" << index1 << ", " << index2
                  << ", " << index3 << "], ";
            ostrm << "native = " << std::fixed << std::showpoint
                  << angle(chain.at(i), chain.at(i+1), chain.at(i+2))
                  << ", ";
            ostrm << "k = " << std::fixed << std::showpoint
                  << this->k_bond_angle_;
            ostrm << "},\n";
        }
        ostrm << "]\n";
    }
    /* dihedral-angle */{
        ostrm << "[[forcefields.local]]\n";
        ostrm << "interaction = \"DihedralAngle\"\n";
        ostrm << "potential   = \"ClementiDihedral\"\n";
        ostrm << "parameters  = [\n";
        for(std::size_t i=0, sz = chain.size() - 3; i<sz; ++i)
        {
            const std::size_t index1 = i + offset;
            const std::size_t index2 = i + offset + 1;
            const std::size_t index3 = i + offset + 2;
            const std::size_t index4 = i + offset + 3;

            ostrm << "{indices = [" << index1 << ", " << index2
                  << ", " << index3 << ", " << index4 << "], ";
            ostrm << "native = " << std::fixed << std::showpoint
                  << dihedral_angle(chain.at(i),   chain.at(i+1),
                                    chain.at(i+2), chain.at(i+3)) << ", ";
            ostrm << "k1 = " << std::fixed << std::showpoint
                  << this->k_dihedral_angle_1_;
            ostrm << ", k3 = " << std::fixed << std::showpoint
                  << this->k_dihedral_angle_3_;
            ostrm << "},\n";
        }
        ostrm << "]\n";
    }
    /* intra-chain-go-contacts */{
        ostrm << "[[forcefields.local]]\n";
        ostrm << "interaction = \"BondLength\"\n";
        ostrm << "potential = \"Go1012Contact\"\n";
        ostrm << "parameters  = [\n";
        for(std::size_t i=0, sz_i = chain.size()-4; i<sz_i; ++i)
        {
            for(std::size_t j=i+4, sz_j = chain.size(); j<sz_j; ++j)
            {
                if(th2 > min_distance_sq_if(
                    chain.at(i)->atoms(), chain.at(j)->atoms(),
                    [](const PDBAtom<coordT>& atom){// ignore hydrogens
                        return atom.atom_name.front() != 'H';
                    }))
                {
                    const std::size_t index1 = i + offset;
                    const std::size_t index2 = j + offset;
                    connections[index1].insert(index2);
                    connections[index2].insert(index1);

                    ostrm << "{indices = [" << index1 << ", " << index2 << "], ";
                    ostrm << "native = " << std::fixed << std::showpoint
                          << distance(chain.at(i), chain.at(j)) << ", ";
                    ostrm << "k = " << std::fixed << std::showpoint
                          << this->k_intra_go_;
                    ostrm << "},\n";
                }
            }
        }
        ostrm << "]\n";
    }
    return connections;
}

template<typename coordT>
bool ClementiGo<coordT>::check_beads_kind(const cg_chain_type& chain) const
{
    bool result = true;
    for(const auto& bead : chain)
    {
        if(bead->kind() != "CarbonAlpha"_str)
        {
            std::cerr << "jarngreipr::ClementiGo::check_beads_kind: Error: \n";
            std::cerr << "    invalid coarse-grained bead kind: "
                      << bead->kind() << '\n';
            std::cerr << "    ClementiGo contains only CarbonAlpha beads.\n";
            result = false;
        }
    }
    return result;
}

}//jarngreipr
#endif /* JARNGREIPR_CLEMENTI_GO */
