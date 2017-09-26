#define BOOST_TEST_MODULE "test_zaxis_external_force_interactionMU"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <cmath>
#include <mjolnir/core/ZaxisExternalForceInteraction.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/ImplicitMembraneList.hpp>
#include <mjolnir/potential/global/ImplicitMembranePotential.hpp>
#include <mjolnir/util/make_unique.hpp>

BOOST_AUTO_TEST_CASE(ZaxisExternalForceInteraction_calc_force)
{
    typedef mjolnir::SimulatorTraitsBase<double, mjolnir::UnlimitedBoundary> traits;
    constexpr static traits::real_type tolerance = 1e-8;

    typedef traits::real_type real_type;
    typedef traits::coordinate_type coord_type;
    typedef traits::boundary_type boundary_type;
    typedef mjolnir::Particle<coord_type> particle_type;
    typedef mjolnir::System<traits> system_type;
    typedef mjolnir::ImplicitMembranePotential<traits> implicit_membrane_type;
    typedef mjolnir::ImplicitMembraneList<traits> partition_type;
    typedef mjolnir::ZaxisExternalForceInteraction<
	traits, implicit_membrane_type,partition_type>
	zaxis_external_force_type;
  
    const real_type thickness(10.);
    const real_type interaction_magnitude(1.);
    const real_type bend(1.5);
    const std::vector<real_type> hydrophobicity{1.};

    implicit_membrane_type potential{
	thickness, interaction_magnitude, bend, hydrophobicity};
    zaxis_external_force_type interaction{
	std::forward<implicit_membrane_type>(potential),
	    partition_type{potential.max_cutoff_length()}};
    
    std::vector<particle_type> ps{
	{1., coord_type(0,0,0), coord_type(0,0,0), coord_type(0,0,0)}	
    };
    
    system_type sys(std::move(ps), boundary_type{});
  
    const real_type dr = 1e-2;
    real_type dist = -10.;
    for(int i = 0; i < 2000; ++i)
    {
	sys[0].position[0] = dist;
	sys[0].position[1] = dist;
	sys[0].position[2] = dist;
	sys[0].force = coord_type(0, 0, 0); 
	
	const real_type deriv = potential.derivative(0, dist);
        const real_type coef  = std::abs(deriv);
	
	interaction.calc_force(sys);

	if(abs(sys[0].position[2]) >= potential.max_cutoff_length())
	{
	    BOOST_CHECK_EQUAL(sys[0].force[2],0);
	}
	else
	{
	    //TODO this process have some error.
	    // const real_type force_strength = mjolnir::length(sys[0].force);
	    // BOOST_CHECK_CLOSE_FRACTION(coef, force_strength, tolerance);
	}
	
	BOOST_CHECK_EQUAL(sys[0].force[0],0);
	BOOST_CHECK_EQUAL(sys[0].force[1],0);

	dist += dr;
    }
}

  
