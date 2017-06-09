#define BOOST_TEST_MODULE "test_global_external_interactionMU"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/GlobalExternalInteractionMU.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/potential/ImplicitMembranePotential.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(GlobalExternal_calc_force)
{
  typedef mjolnir::SimulatorTraitsBase<double, mjolnir::Unlimited> traits;
  constexpr static traits::real_type tolerance = 1e-8;

  typedef traits::real_type real_type;
  typedef traits::coordinate_type coord_type;
  typedef traits::boundary_type boudary_type;
  typedef mjolnir::Particle<coord_type> particle_type;
  typedef mjolnir::System<traits> system_type;
  typedef mjolnir::ImplicitMembranePotential<traits> implicit_membrane_type;
  typedef mjolnir::GlobalExternalInteraction<traits, implicit_membrane_type>
    global_external_type;

  std::vector<particle_type> ps{
    {1., coord_type(0,0,0), coord_type(0,0,0), coord_type(0,0,0)}
  };
  system_type sys(std::move(ps), boundary_type{});

  
  
  implicit_membrane_type potential{  };
  
  
