#define BOOST_TEST_MODULE "test_dihedral_angle_utility_function"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/interaction/local/DihedralAngleInteraction.hpp>
#include <mjolnir/interaction/local/dihedral_angle.hpp>
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

BOOST_AUTO_TEST_CASE(DihedralAngleInteraction_and_calc_dihedral_angle_force)
{
    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using potential_type   = mjolnir::ConstantPotential<real_type>;
    using interaction_type = mjolnir::DihedralAngleInteraction<traits_type, potential_type>;

    const auto tol = boost::test_tools::tolerance(1e-5);

    system_type sys(4, boundary_type{});

    // -----------------------------------------------------------------------
    // check calc_angle calculate the correct angle

    const int N = 1800;
    const real_type dtheta = mjolnir::math::constants<real_type>::pi / N;
    for(int i = -1800; i < 1800; ++i)
    {
        const auto theta = i * dtheta;

        const auto p0 = mjolnir::math::make_coordinate<coord_type>(1, 0, 1);
        const auto p1 = mjolnir::math::make_coordinate<coord_type>(0, 0, 1);
        const auto p2 = mjolnir::math::make_coordinate<coord_type>(0, 0, 0);
        const auto p3 = mjolnir::math::make_coordinate<coord_type>(
                std::cos(theta), -std::sin(theta), 0.0);

        const auto angle = mjolnir::calc_dihedral_angle(sys, p0, p1, p2, p3);
        BOOST_TEST(theta == angle, tol);

        const auto f = mjolnir::calc_dihedral_angle_force(sys, p0, p1, p2, p3);

        BOOST_TEST(theta == f.first, tol);
        BOOST_TEST(angle == f.first, tol);

        const auto dot_v01_f0 = mjolnir::math::dot_product(f.second.at(0), p0 - p1);
        const auto dot_v12_f0 = mjolnir::math::dot_product(f.second.at(0), p1 - p2);

        const auto dot_v23_f3 = mjolnir::math::dot_product(f.second.at(3), p2 - p3);
        const auto dot_v12_f3 = mjolnir::math::dot_product(f.second.at(3), p1 - p2);

        BOOST_TEST(dot_v01_f0 == 0.0, tol);
        BOOST_TEST(dot_v12_f0 == 0.0, tol);
        BOOST_TEST(dot_v12_f3 == 0.0, tol);
        BOOST_TEST(dot_v23_f3 == 0.0, tol);

        const auto f_sum = f.second.at(0) + f.second.at(1) + f.second.at(2) + f.second.at(3);

        // sum of the forces should be zero, no net translation or rotation
        BOOST_TEST(mjolnir::math::X(f_sum) == 0.0, tol);
        BOOST_TEST(mjolnir::math::Y(f_sum) == 0.0, tol);
        BOOST_TEST(mjolnir::math::Z(f_sum) == 0.0, tol);
    }

    // -----------------------------------------------------------------------

    potential_type   potential;
    interaction_type interaction("none", {{ {{0, 1, 2, 3}}, potential}});

    sys.at(0).mass = 1.0;
    sys.at(1).mass = 1.0;
    sys.at(2).mass = 1.0;
    sys.at(3).mass = 1.0;
    sys.at(0).rmass = 1.0;
    sys.at(1).rmass = 1.0;
    sys.at(2).rmass = 1.0;
    sys.at(3).rmass = 1.0;

    sys.at(0).position = coord_type(1.0, 0.0, 0.0);
    sys.at(1).position = coord_type(0.0, 0.0, 1.0);
    sys.at(2).position = coord_type(1.0, 1.0, 1.0);
    sys.at(0).velocity = coord_type(0.0, 0.0, 0.0);
    sys.at(1).velocity = coord_type(0.0, 0.0, 0.0);
    sys.at(2).velocity = coord_type(0.0, 0.0, 0.0);
    sys.at(3).velocity = coord_type(0.0, 0.0, 0.0);
    sys.at(0).force    = coord_type(0.0, 0.0, 0.0);
    sys.at(1).force    = coord_type(0.0, 0.0, 0.0);
    sys.at(2).force    = coord_type(0.0, 0.0, 0.0);
    sys.at(3).force    = coord_type(0.0, 0.0, 0.0);

    sys.at(0).name  = "X";
    sys.at(1).name  = "X";
    sys.at(2).name  = "X";
    sys.at(3).name  = "X";
    sys.at(0).group = "NONE";
    sys.at(1).group = "NONE";
    sys.at(2).group = "NONE";
    sys.at(3).group = "NONE";

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(0.0, 1.0);

    for(int i = 1; i < N; ++i)
    {
        sys[0].force = mjolnir::math::make_coordinate<coord_type>(0, 0, 0);
        sys[1].force = mjolnir::math::make_coordinate<coord_type>(0, 0, 0);
        sys[2].force = mjolnir::math::make_coordinate<coord_type>(0, 0, 0);
        sys[3].force = mjolnir::math::make_coordinate<coord_type>(0, 0, 0);

        sys[0].position = mjolnir::math::make_coordinate<coord_type>(-1.5, 0, 0);
        sys[1].position = mjolnir::math::make_coordinate<coord_type>(-uni(mt), -uni(mt), -uni(mt));
        sys[2].position = mjolnir::math::make_coordinate<coord_type>( uni(mt),  uni(mt),  uni(mt));
        sys[3].position = mjolnir::math::make_coordinate<coord_type>(   0, 0, 1.5);

        const auto a = mjolnir::calc_dihedral_angle      (sys, sys.position(0), sys.position(1), sys.position(2), sys.position(3));
        const auto f = mjolnir::calc_dihedral_angle_force(sys, sys.position(0), sys.position(1), sys.position(2), sys.position(3));

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

        BOOST_TEST(mjolnir::math::X(sys[3].force) == mjolnir::math::X(f.second[3]), tol);
        BOOST_TEST(mjolnir::math::Y(sys[3].force) == mjolnir::math::Y(f.second[3]), tol);
        BOOST_TEST(mjolnir::math::Z(sys[3].force) == mjolnir::math::Z(f.second[3]), tol);
    }
}
