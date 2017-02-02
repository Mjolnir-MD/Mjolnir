#ifndef JARNGREIPR_CLEMENTI_GO
#define JARNGREIPR_CLEMENTI_GO
#include "Model.hpp"
#include "CarbonAlpha.hpp"
#include <jarngreipr/geometry/distance.hpp>
#include <jarngreipr/geometry/angle.hpp>
#include <jarngreipr/geometry/dihedral.hpp>
#include <iterator>
#include <vector>

namespace jarngreipr
{

template<typename traitsT>
class ClementiGo : public Model<traitsT>
{
  public:

    typedef Model<traitsT> base_type;
    typedef typename base_type::real_type            real_type;
    typedef typename base_type::coordinate_type      coordinate_type;
    typedef typename base_type::bead_type            bead_type;
    typedef typename base_type::bead_container_type  bead_container_type;
    typedef typename base_type::chain_type           chain_type;
    typedef typename base_type::chain_container_type chain_container_type;
    typedef std::vector<std::size_t>                 index_list_type;
    typedef std::vector<index_list_type>             exception_list_type;
    template<std::size_t N>
    using interaction_type = typename base_type::interaction_type<N>;

  public:

    ClementiGo(): contact_threshold_(6.5){}
    explicit ClementiGo(const real_type th): contact_threshold_(th){}
    ~ClementiGo() = default;

    void make(const chain_container_type& chain) override;

    real_type& contact_threshold()       {return contact_threshold_;}
    real_type  contact_threshold() const {return contact_threshold_;}

    interaction_type<2>&       bond()            {return bond_length_;}
    interaction_type<2> const& bond()      const {return bond_length_;}
    interaction_type<2>&       contact()         {return go_contact_;}
    interaction_type<2> const& contact()   const {return go_contact_;}
    interaction_type<3>&       angle()           {return bond_angle_;}
    interaction_type<3> const& angle()     const {return bond_angle_;}
    interaction_type<4>&       dihedral()        {return dihedral_angle_;}
    interaction_type<4> const& dihedral()  const {return dihedral_angle_;}
    exception_list_type &      exception()       {return exception_;}
    exception_list_type const& exception() const {return exception_;}

  private:

    real_type contact_threshold_;

    std::vector<interaction_type<2>> bond_length_;
    std::vector<interaction_type<2>> go_contact_;
    std::vector<interaction_type<3>> bond_angle_;
    std::vector<interaction_type<4>> dihedral_angle_;
    exception_list_type exception_; // !< for EXV exception
};

template<typename traitsT>
void ClementiGo<traitsT>::make(const chain_container_type& chains)
{
    const real_type threshold2 = contact_threshold_ * contact_threshold_;

    for(auto chain = chains.cbegin(); chain != chains.cend(); ++chain)
    {
        const std::size_t cg_begin = this->beads_.size();
        // make cg-beads
        for(auto iter = chain->cbegin(); iter != chain->cend(); ++iter)
            this->beads_.emplace_back(*iter);

        { // ditect go-contact
        std::size_t i = cg_begin;
        exception_.resize(this->beads.size());
        for(auto iter = chain->cbegin(); iter != chain->cend()-3; ++iter)
        {
             std::size_t j = i+3;
             for(auto jter = iter + 3; jter != chain->cend(); ++jter)
             {
                if(min_distance_sq(*iter, *jter) < threshold2)
                {
                    const auto dist = distance(this->beads_.at(i).position(),
                                               this->beads_.at(j).position());
                    this->go_contact_.emplace_back({{i, j}}, dist);
                    this->exception_.at(i).push_back(j);
                    this->exception_.at(j).push_back(i);
                }
                ++j;
             }
             ++i;
        }
        } // go-contact

        {// bond, angle, dihd
        for(std::size_t i=cg_begin; i<beads_.size(); ++i)
        {
            this->exception_.at(i).reserve(8);
            this->exception_.at(i).push_back(i);

            if(i+1 < beads_.size()) // bond
            {
                bond_length_.emplace_back({{i, i+1}}, distance(
                        this->beads_.at(i).position(),
                        this->beads_.at(i+1).position()));
                for(std::size_t j=i; j<=i+1; ++j)
                for(std::size_t k=i; k<=i+1; ++k)
                    if(j != k) this->exception_.at(j).push_back(k);
            }

            if(i+2 < beads_.size()) // angle
            {
                bond_angle_.emplace_back({{i, i+1, i+2}}, angle(
                        this->beads_.at(i).position(),
                        this->beads_.at(i+1).position(),
                        this->beads_.at(i+2).position()));
                for(std::size_t j=i; j<=i+2; ++j)
                for(std::size_t k=i; k<=i+2; ++k)
                    if(j != k) this->exception_.at(j).push_back(k);
            }

            if(i+3 < beads_.size()) // dihd
            {
                dihedral_angle_.emplace_back({{i, i+1, i+2, i+3}},
                        dihedral(this->beads_.at(i  ).position(),
                                 this->beads_.at(i+1).position(),
                                 this->beads_.at(i+2).position(),
                                 this->beads_.at(i+3).position()));
                for(std::size_t j=i; j<=i+3; ++j)
                for(std::size_t k=i; k<=i+3; ++k)
                    if(j != k) this->exception_.at(j).push_back(k);
            }
        }
        }// bond, angle, dihd
    }

    // inter-chain go-contact assuming(num of residue == num of Ca)
    std::size_t i_begin = 0;
    std::size_t j_begin = 0;
    for(auto iter = chains.cbegin(); iter != chains.cend()-1; ++iter)
    {
        j_begin = i_begin + iter->size();
        for(auto jter = iter+1; jter != chains.cend(); ++jter)
        {
            std::size_t i = i_begin;
            std::size_t j = j_begin;
            for(auto lhs = iter->cbegin(); lhs != iter->cend(); ++lhs)
            {// PDBResidue
                for(auto rhs = iter->cbegin(); rhs != iter->cend(); ++rhs)
                {
                    if(min_distance_sq(*lhs, *rhs) < threshold2)
                    {
                        const auto dist = distance(this->beads_.at(i).position(),
                                                   this->beads_.at(j).position());
                        this->go_contact_.emplace_back({{i, j}}, dist);
                        this->exception_.at(i).push_back(j);
                        this->exception_.at(j).push_back(i);
                    }
                    ++j;
                }
                ++i;
            }
            j_begin += jter->size();
        }
        i_begin += iter->size();
    }


    return;
}

}//jarngreipr
#endif /* JARNGREIPR_CLEMENTI_GO */
