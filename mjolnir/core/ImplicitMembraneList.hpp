#ifndef MJOLNIR_IMPLICIT_MEMBRANE_LIST
#define MJOLNIR_IMPLICIT_MEMBRANE_LIST

#include <vector>
#include <cmath>

namespace mjolnir
{
template<typename traitT>
class ImplicitMembraneList
{
  public:
    typedef traitT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef typename std::vector<std::size_t> obj_index_array;
    
    ImplicitMembraneList() = default;
    ~ImplicitMembraneList() = default;

    ImplicitMembraneList(const real_type cutoff) noexcept :cutoff_(cutoff) {}
    
    void initialize(const system_type& sys);
    void update(const system_type& sys);
    obj_index_array const& obj_indices()const noexcept {return obj_indices_;}
    
  private:

    real_type cutoff_;
    obj_index_array obj_indices_;
};

template<typename traitT>
void ImplicitMembraneList<traitT>::initialize(const system_type& sys)
{
    for(std::size_t particle_idx = 0;  particle_idx < sys.size(); ++particle_idx)
    {
	if(std::abs(sys[particle_idx].position[2]) < cutoff_)
	    obj_indices_.push_back(particle_idx);
    }
    return;
}

template<typename traitT>
void ImplicitMembraneList<traitT>::update(const system_type& sys)
{
    obj_indices_.clear();
    initialize(sys);
    return;
}

}

#endif /* MJOLNIR_IMPLICIT_MEMBRANE_LIST */
