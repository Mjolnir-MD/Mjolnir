#define BOOST_TEST_MODULE "test_bond_angle_utility_function"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/interaction/local/BondAngleInteraction.hpp>
#include <mjolnir/interaction/local/bond_angle.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <random>

namespace mjolnir
{
template<typename realT>
class ConstantPotential
{
  public:
    using real_type = realT;

  public:
    ConstantPotential()  = default;
    ~ConstantPotential() = default;

    real_type potential(const real_type v) const noexcept
    {
        return -v;
    }

    real_type derivative(const real_type) const noexcept
    {
        return real_type(-1.0);
    }

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() noexcept {return "Constant";}

    real_type cutoff() const noexcept // no cutoff exists.
    {
        return std::numeric_limits<real_type>::infinity();
    }
};
} // mjolnir

BOOST_AUTO_TEST_CASE(BondAngleInteraction_and_calc_bond_angle_force)
{
    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::ConstantPotential<real_type>;
    using interaction_type = mjolnir::BondAngleInteraction<traits_type, potential_type>;

    const auto tol = boost::test_tools::tolerance(1e-5);

    system_type sys(3, boundary_type{});

    // -----------------------------------------------------------------------
    // check calc_angle calculate the correct angle

    const int N = 1800;
    const real_type dtheta = mjolnir::math::constants<real_type>::pi / N;
    for(int i = 1; i < N; ++i)
    {
        const auto theta = i * dtheta;

        const auto p0 = mjolnir::math::make_coordinate<coord_type>(1, 0, 0);
        const auto p1 = mjolnir::math::make_coordinate<coord_type>(0, 0, 0);
        const auto p2 = mjolnir::math::make_coordinate<coord_type>(
                std::cos(theta), std::sin(theta), 0);

        const auto angle = mjolnir::calc_bond_angle(sys, p0, p1, p2);
        BOOST_TEST(theta == angle, tol);

        const auto f = mjolnir::calc_bond_angle_force(sys, p0, p1, p2);

        BOOST_TEST(theta == f.first, tol);
        BOOST_TEST(angle == f.first, tol);

        const auto dot_v10_f = mjolnir::math::dot_product(f.second.at(0), p0 - p1);
        const auto dot_v12_f = mjolnir::math::dot_product(f.second.at(2), p2 - p1);

        const auto f_sum = f.second.at(0) + f.second.at(1) + f.second.at(2);

        BOOST_TEST(dot_v10_f           == 0.0, tol);
        BOOST_TEST(dot_v12_f           == 0.0, tol);

        // each force parpendicular to Z axis
        BOOST_TEST(mjolnir::math::Z(f.second.at(0)) == 0.0, tol);
        BOOST_TEST(mjolnir::math::Z(f.second.at(1)) == 0.0, tol);
        BOOST_TEST(mjolnir::math::Z(f.second.at(2)) == 0.0, tol);

        // sum of the forces should be zero, no net translation or rotation
        BOOST_TEST(mjolnir::math::X(f_sum) == 0.0, tol);
        BOOST_TEST(mjolnir::math::Y(f_sum) == 0.0, tol);
        BOOST_TEST(mjolnir::math::Z(f_sum) == 0.0, tol);
    }

    // -----------------------------------------------------------------------

    potential_type   potential;
    interaction_type interaction("none", {{ {{0,1,2}}, potential}});

    sys.at(0).mass = 1.0;
    sys.at(1).mass = 1.0;
    sys.at(2).mass = 1.0;
    sys.at(0).rmass = 1.0;
    sys.at(1).rmass = 1.0;
    sys.at(2).rmass = 1.0;

    sys.at(0).position = coord_type(1.0, 0.0, 0.0);
    sys.at(1).position = coord_type(0.0, 0.0, 1.0);
    sys.at(2).position = coord_type(1.0, 1.0, 1.0);
    sys.at(0).velocity = coord_type(0.0, 0.0, 0.0);
    sys.at(1).velocity = coord_type(0.0, 0.0, 0.0);
    sys.at(2).velocity = coord_type(0.0, 0.0, 0.0);
    sys.at(0).force    = coord_type(0.0, 0.0, 0.0);
    sys.at(1).force    = coord_type(0.0, 0.0, 0.0);
    sys.at(2).force    = coord_type(0.0, 0.0, 0.0);

    sys.at(0).name  = "X";
    sys.at(1).name  = "X";
    sys.at(2).name  = "X";
    sys.at(0).group = "NONE";
    sys.at(1).group = "NONE";
    sys.at(2).group = "NONE";

    const auto init = sys;

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(0.0, 1.0);

    for(int i = 1; i < N; ++i)
    {
        sys[0].force = mjolnir::math::make_coordinate<coord_type>(0, 0, 0);
        sys[1].force = mjolnir::math::make_coordinate<coord_type>(0, 0, 0);
        sys[2].force = mjolnir::math::make_coordinate<coord_type>(0, 0, 0);

        sys[0].position = mjolnir::math::make_coordinate<coord_type>(-0.5, -0.5, -0.5);
        sys[1].position = mjolnir::math::make_coordinate<coord_type>(uni(mt), uni(mt), uni(mt));
        sys[2].position = mjolnir::math::make_coordinate<coord_type>( 1.5,  1.5,  1.5);

        const auto a = mjolnir::calc_bond_angle      (sys, sys.position(0), sys.position(1), sys.position(2));
        const auto f = mjolnir::calc_bond_angle_force(sys, sys.position(0), sys.position(1), sys.position(2));

        BOOST_TEST(a == f.first, tol);

        interaction.calc_force(sys);

        // perpendicular to z axis
        BOOST_TEST(mjolnir::math::X(sys[0].force) == mjolnir::math::X(f.second[0]), tol);
        BOOST_TEST(mjolnir::math::Y(sys[0].force) == mjolnir::math::Y(f.second[0]), tol);
        BOOST_TEST(mjolnir::math::Z(sys[0].force) == mjolnir::math::Z(f.second[0]), tol);

        BOOST_TEST(mjolnir::math::X(sys[1].force) == mjolnir::math::X(f.second[1]), tol);
        BOOST_TEST(mjolnir::math::Y(sys[1].force) == mjolnir::math::Y(f.second[1]), tol);
        BOOST_TEST(mjolnir::math::Z(sys[1].force) == mjolnir::math::Z(f.second[1]), tol);

        BOOST_TEST(mjolnir::math::X(sys[2].force) == mjolnir::math::X(f.second[2]), tol);
        BOOST_TEST(mjolnir::math::Y(sys[2].force) == mjolnir::math::Y(f.second[2]), tol);
        BOOST_TEST(mjolnir::math::Z(sys[2].force) == mjolnir::math::Z(f.second[2]), tol);
    }
}
