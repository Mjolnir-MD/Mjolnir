#ifndef JARNGREIPR_CLEMENTI_GO
#define JARNGREIPR_CLEMENTI_GO
#include "Model.hpp"
#include "CarbonAlpha.hpp"
#include <jarngreipr/geometry/distance.hpp>
#include <jarngreipr/geometry/angle.hpp>
#include <jarngreipr/geometry/dihedral.hpp>
#include <mjolnir/util/make_unique.hpp>
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
    using interaction_type = std::pair<std::array<std::size_t, N>, real_type>;

  public:

    ClementiGo(): contact_threshold_(6.5){}
    explicit ClementiGo(const real_type th): contact_threshold_(th){}
    ~ClementiGo() = default;

    void make(const chain_container_type& chain) override;

    real_type& contact_threshold()       {return contact_threshold_;}
    real_type  contact_threshold() const {return contact_threshold_;}

    std::vector<interaction_type<2>>&       bonds()            {return bond_length_;}
    std::vector<interaction_type<2>> const& bonds()      const {return bond_length_;}
    std::vector<interaction_type<2>>&       contacts()         {return go_contact_;}
    std::vector<interaction_type<2>> const& contacts()   const {return go_contact_;}
    std::vector<interaction_type<3>>&       angles()           {return bond_angle_;}
    std::vector<interaction_type<3>> const& angles()     const {return bond_angle_;}
    std::vector<interaction_type<4>>&       dihedrals()        {return dihedral_angle_;}
    std::vector<interaction_type<4>> const& dihedrals()  const {return dihedral_angle_;}
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
            this->beads_.emplace_back(mjolnir::make_unique<CarbonAlpha<traitsT>>(
                        iter->atoms(), iter->residue_name()));

        { // ditect go-contact
        std::size_t i = cg_begin;
        this->exception_.resize(this->beads_.size());
        for(auto iter = chain->cbegin(); iter != chain->cend()-4; ++iter)
        {
             std::size_t j = i+4;
             for(auto jter = iter + 4; jter != chain->cend(); ++jter)
             {
                if(min_distance_sq(*iter, *jter) < threshold2)
                {
                    const auto dist = distance(this->beads_.at(i)->position(0),
                                               this->beads_.at(j)->position(0));
                    std::array<std::size_t, 2> indices{{i, j}};
                    this->go_contact_.emplace_back(indices, dist);
                    this->exception_.at(i).push_back(j);
                    this->exception_.at(j).push_back(i);
                }
                ++j;
             }
             ++i;
        }
        } // go-contact

        {// bond, angle, dihd
        for(std::size_t i=cg_begin; i<this->beads_.size(); ++i)
        {
            this->exception_.at(i).reserve(8);
            for(std::size_t j = std::max(static_cast<int>(i)-3, 0);
                    j <= std::min(i+3, this->beads_.size()); ++j)
            {
                this->exception_.at(i).push_back(j);
            }

            if(i+1 < this->beads_.size()) // bond
            {
                const std::array<std::size_t, 2> indices{{i, i+1}};
                bond_length_.emplace_back(indices, distance(
                        this->beads_.at(i)->position(0),
                        this->beads_.at(i+1)->position(0)));
            }

            if(i+2 < this->beads_.size()) // angle
            {
                const std::array<std::size_t, 3> indices{{i, i+1, i+2}};
                bond_angle_.emplace_back(indices, angle(
                        this->beads_.at(i)->position(0),
                        this->beads_.at(i+1)->position(0),
                        this->beads_.at(i+2)->position(0)));
            }

            if(i+3 < this->beads_.size()) // dihd
            {
                const std::array<std::size_t, 4> indices{{i, i+1, i+2, i+3}};
                dihedral_angle_.emplace_back(indices, dihedral_angle(
                            this->beads_.at(i  )->position(0),
                            this->beads_.at(i+1)->position(0),
                            this->beads_.at(i+2)->position(0),
                            this->beads_.at(i+3)->position(0)));
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
            for(auto lhs = iter->cbegin(); lhs != iter->cend(); ++lhs)
            {
                std::size_t j = j_begin;
                for(auto rhs = jter->cbegin(); rhs != jter->cend(); ++rhs)
                {
                    if(min_distance_sq(*lhs, *rhs) < threshold2)
                    {
                        const auto dist = distance(this->beads_.at(i)->position(0),
                                                   this->beads_.at(j)->position(0));
                        const std::array<std::size_t, 2> indices{{i, j}};
                        this->go_contact_.emplace_back(indices, dist);
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
