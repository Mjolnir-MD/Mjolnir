#ifndef JARNGREIPR_GO_CONTACT
#define JARNGREIPR_GO_CONTACT
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
class GoContact final : public InterChainForceFieldGenerator<coordT>
{
  public:
    typedef coordT coordinate_type;
    typedef InterChainForceFieldGenerator<coordT> base_type;
    typedef typename base_type::atom_type         atom_type;
    typedef typename base_type::residue_type      residue_type;
    typedef typename base_type::chain_type        chain_type;
    typedef typename base_type::cg_chain_type     cg_chain_type;
    typedef typename base_type::connection_info   connection_info;
    typedef typename scalar_type_of<coordinate_type>::type real_type;

  public:

    GoContact()
        : contact_threshold_(6.5), k_(0.3)
    {}
    ~GoContact() override = default;

    explicit GoContact(const real_type contact_threshold)
        : contact_threshold_(contact_threshold), k_(0.3)
    {}

    GoContact(const real_type contact_threshold, const real_type k)
        : contact_threshold_(contact_threshold), k_(k)
    {}

    connection_info
    generate(std::ostream& ostrm,
        const std::vector<std::pair<cg_chain_type, std::size_t>>& chain_offset
        ) const;

    bool check_beads_kind(const cg_chain_type& chain) const override
    {
        return true; // accept all kinds of beads.
    }

    real_type& contact_threshold()       noexcept {return contact_threshold_;}
    real_type  contact_threshold() const noexcept {return contact_threshold_;}
    real_type& k()       noexcept {return k_;}
    real_type  k() const noexcept {return k_;}

  private:
    real_type contact_threshold_;
    real_type k_;
};

template<typename coordT>
typename GoContact<coordT>::connection_info
GoContact<coordT>::generate(std::ostream& ostrm,
        const std::vector<std::pair<cg_chain_type, std::size_t>>& chain_offset
        ) const
{
    for(const auto& ch_ofs : chain_offset)
    {
        if(false == this->check_bead_kind(ch_ofs.first))
        {
            // never reach here...; for consistency.
            throw std::invalid_argument(
                    "jarngreipr::GoContact::generate: "_str +
                    "invalid bead kind appear in argument `chain`."_str);
        }
    }

    connection_info connections;
    for(const auto& ch_ofs : chain_offset)
    {
        const auto& chain  = ch_ofs.first;
        const auto  offset = ch_ofs.second;
        for(std::size_t i=0; i < chain.first.size(); ++i)
        {
            const std::size_t index = i + offset;
            connections[index].insert(index);
        }
    }

    ostrm << "[[forcefields.local]]\n";
    ostrm << "interaction = \"BondLength\"\n";
    ostrm << "potential = \"Go1012Contact\"\n";
    ostrm << "parameters  = [\n";

    const real_type th2 = this->contact_threshold_ * this->contact_threshold_;
    for(std::size_t chain1_index = 0, chain1_end = chain_offset.size()-1;
            chain1_index < chain1_end; ++chain1_index)
    {
        const auto& chain1        = chain_offset.at(chain1_index).first;
        const auto& chain1_offset = chain_offset.at(chain1_index).second;

        for(std::size_t chain2_index = chain1_index+1;
                chain2_index < chain_offset.size(); ++chain2_index)
        {
            const auto& chain2        = chain_offset.at(chain2_index).first;
            const auto& chain2_offset = chain_offset.at(chain2_index).second;

            for(std::size_t i=0, i < chain1.size(); ++i)
            {
                for(std::size_t j=0, j < chain2.size(); ++j)
                {
                    if(th2 > min_distance_sq_if(
                        chain1.at(i)->atoms(), chain2.at(j).atoms(),
                        [](const PDBAtom<coordT>& atom){// ignore hydrogens
                            return atom.atom_name.front() != 'H';
                        }))
                    {
                        connections[i].insert(j);
                        connections[j].insert(i);

                        ostrm << "{indices = [" << i << ", " << j << "], ";
                        ostrm << "native = " << std::fixed << std::showpoint
                              << distance(chain1.at(i), chain2.at(j)) << ", ";
                        ostrm << "k = " << std::fixed << std::showpoint
                              << this->k_;
                        ostrm << "},\n";
                    }
                }
            }
        }
    }
    ostrm << "]\n";
    return connections;
}

} // mjolnir
#endif// JARNGREIPR_GO_CONTACT
