#ifndef JARNGREIPR_CLEMENTI_GO
#define JARNGREIPR_CLEMENTI_GO
#include "Model.hpp"
#include "CarbonAlpha.hpp"
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

    ClementiGo() = default;
    ClementiGo(const chain_container_type& chain)
    {
        this->make(chain);
    }
    ~ClementiGo() = default;

    void make(const chain_container_type& chain) override;

    interaction_type<2>&       bond()       {return bond_length_;}
    interaction_type<2> const& bond() const {return bond_length_;}

    interaction_type<2>&       contact()       {return go_contact_;}
    interaction_type<2> const& contact() const {return go_contact_;}

    interaction_type<3>&       angle()       {return bond_angle_;}
    interaction_type<3> const& angle() const {return bond_angle_;}

    interaction_type<4>&       dihedral()       {return dihedral_angle_;}
    interaction_type<4> const& dihedral() const {return dihedral_angle_;}

    exception_list_type &      exception()       {return exception_;}
    exception_list_type const& exception() const {return exception_;}

  private:

    interaction_type<2> bond_length_;
    interaction_type<2> go_contact_;
    interaction_type<3> bond_angle_;
    interaction_type<4> dihedral_angle_;
    exception_list_type exception_; // !< for EXV exception
};

template<typename traitsT>
void ClementiGo<traitsT>::make(const chain_container_type& chain)
{
    ;// TODO
    return;
}




}//jarngreipr
#endif /* JARNGREIPR_CLEMENTI_GO */
