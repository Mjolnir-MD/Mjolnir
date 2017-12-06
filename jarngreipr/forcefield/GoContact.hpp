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
             const std::vector<cg_chain_type>& chains) const;

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
        const std::vector<cg_chain_type>& chains) const
{
    for(const auto& chain : chains)
    {
        if(false == this->check_bead_kind(chain))
        {
            // never reach here...; for consistency.
            throw std::invalid_argument(
                    "jarngreipr::GoContact::generate: "_str +
                    "invalid bead kind appear in argument `chain`."_str);
        }
    }

    connection_info connections;
    for(const auto& chain : chains)
    {
        for(std::size_t i=0; i < chain.size(); ++i)
        {
            const auto& b = chain.at(i)
            const std::size_t index = b->index();
            connections[index].insert(index);
        }
    }

    ostrm << "[[forcefields.local]]\n";
    ostrm << "interaction = \"BondLength\"\n";
    ostrm << "potential = \"Go1012Contact\"\n";
    ostrm << "parameters  = [\n";

    const real_type th2 = this->contact_threshold_ * this->contact_threshold_;
    for(std::size_t ch1idx = 0, ch1end = chains.size()-1;
            ch1idx < ch1end; ++ch1idx)
    {
        const auto& chain1 = chains.at(ch1idx);

        for(std::size_t ch2idx = ch1idx+1, ch2end = chains.size();
                ch2idx < ch2end; ++ch2idx)
        {
            const auto& chain2 = chains.at(ch2idx);

            for(std::size_t i=0, i < chain1.size(); ++i)
            {
                for(std::size_t j=0, j < chain2.size(); ++j)
                {
                    const auto& bead1 = chain1.at(i);
                    const auto& bead2 = chain2.at(j);
                    if(th2 > min_distance_sq_if(bead1->atoms(), bead2->atoms(),
                        [](const PDBAtom<coordT>& atom){// ignore hydrogens
                            return atom.atom_name.front() != 'H';
                        }))
                    {
                        const std::size_t idx1 = bead1->index();
                        const std::size_t idx2 = bead2->index();
                        connections[idx1].insert(idx2);
                        connections[idx2].insert(idx1);

                        ostrm << "{indices = [" << idx1 << ", " << idx2 << "], ";
                        ostrm << "native = " << std::fixed << std::showpoint
                              << distance(bead1, bead2) << ", ";
                        ostrm << "k = " << std::fixed << std::showpoint << this->k_;
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
