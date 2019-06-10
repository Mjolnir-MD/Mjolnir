#define BOOST_TEST_MODULE "test_bond_angle_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/interaction/local/BondAngleInteraction.hpp>
#include <mjolnir/potential/local/HarmonicPotential.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/util/make_unique.hpp>

#include <random>

BOOST_AUTO_TEST_CASE(BondAngleInteraction_force)
{
    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type       = traits_type::real_type;
    using coord_type      = traits_type::coordinate_type;
    using boundary_type   = traits_type::boundary_type;
    using system_type     = mjolnir::System<traits_type>;
    using harmonic_type   = mjolnir::HarmonicPotential<real_type>;
    using bond_angle_type = mjolnir::BondAngleInteraction<traits_type, harmonic_type>;

    constexpr real_type tol = 1e-7;

    auto normalize = [](const coord_type& v){return v / mjolnir::math::length(v);};

    const real_type k(1e0);
    const real_type native(mjolnir::math::constants<real_type>::pi * 2.0 / 3.0); // 120 degree

    harmonic_type potential{k, native};
    bond_angle_type interaction("none", {{ {{0,1,2}}, potential}});

    const coord_type pos1(1., 0., 0.);
    const coord_type pos2(0., 0., 0.);
    system_type sys(3, boundary_type{});

    sys.at(0).mass = 1.0;
    sys.at(1).mass = 1.0;
    sys.at(2).mass = 1.0;
    sys.at(0).rmass = 1.0;
    sys.at(1).rmass = 1.0;
    sys.at(2).rmass = 1.0;

    sys.at(0).position = pos1;
    sys.at(1).position = pos2;
    sys.at(2).position = coord_type(0,0,0);
    sys.at(0).velocity = coord_type(0,0,0);
    sys.at(1).velocity = coord_type(0,0,0);
    sys.at(2).velocity = coord_type(0,0,0);
    sys.at(0).force    = coord_type(0,0,0);
    sys.at(1).force    = coord_type(0,0,0);
    sys.at(2).force    = coord_type(0,0,0);

    sys.at(0).name  = "X";
    sys.at(1).name  = "X";
    sys.at(2).name  = "X";
    sys.at(0).group = "NONE";
    sys.at(1).group = "NONE";
    sys.at(2).group = "NONE";

    const int N = 1800;
    const real_type dtheta = mjolnir::math::constants<real_type>::pi  / N;
    for(int i = 1; i < N; ++i)
    {
        BOOST_TEST(mjolnir::math::length(sys[0].position - pos1) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[1].position - pos2) == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[0].velocity)        == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[1].velocity)        == 0.0, boost::test_tools::tolerance(tol));
        sys[0].force = coord_type(0,0,0);
        sys[1].force = coord_type(0,0,0);
        sys[2].force = coord_type(0,0,0);

        const real_type theta = i * dtheta;
        const coord_type pos3(std::cos(theta), std::sin(theta), 0e0);
        sys[2].position = pos3;

        const real_type deriv = potential.derivative(theta);
        const real_type coef = std::abs(deriv);

        interaction.calc_force(sys);

        // magnitude
        // if radius == 1e0, then force strength is equal to dV/dtheta.
        BOOST_TEST(mjolnir::math::length(sys[0].position - sys[1].position) == 1.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::length(sys[2].position - sys[1].position) == 1.0, boost::test_tools::tolerance(tol));

        const real_type force_strength1 = mjolnir::math::length(sys[0].force);
        const real_type force_strength3 = mjolnir::math::length(sys[2].force);
        if(i == 1200) // most stable point
        {
            BOOST_TEST(coef == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(coef == 0.0, boost::test_tools::tolerance(tol));
        }
        else
        {
            BOOST_TEST(coef == force_strength1, boost::test_tools::tolerance(tol));
            BOOST_TEST(coef == force_strength3, boost::test_tools::tolerance(tol));
        }

        // force applied to center particle is equal to sum of others
        const coord_type sum = sys[0].force + sys[1].force + sys[2].force;
        BOOST_TEST(mjolnir::math::length(sum) == 0.0, boost::test_tools::tolerance(tol));

        // direction
        if(i == 1200) // most stable point
        {
            BOOST_TEST(std::abs(force_strength1) == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(std::abs(force_strength3) == 0.0, boost::test_tools::tolerance(tol));
        }
        else if(i < 1200) // narrow
        {
            // perpendicular to radius vector
            const real_type dot1 = mjolnir::math::dot_product(sys[0].force, sys[0].position - sys[1].position);
            const real_type dot2 = mjolnir::math::dot_product(sys[2].force, sys[2].position - sys[1].position);
            BOOST_TEST(dot1 == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dot2 == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f1(0., -1., 0.);
            BOOST_TEST(mjolnir::math::length(normalize(sys[0].force) - f1) == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f3(
                    cos(theta + mjolnir::math::constants<real_type>::pi / 2.),
                    sin(theta + mjolnir::math::constants<real_type>::pi / 2.),
                    0.);
            BOOST_TEST(mjolnir::math::length(normalize(sys[2].force) - f3) == 0.0, boost::test_tools::tolerance(tol));
        }
        else if(i > 1200) // extensive
        {
            const real_type dot1 = mjolnir::math::dot_product(sys[0].force, sys[0].position - sys[1].position);
            const real_type dot2 = mjolnir::math::dot_product(sys[2].force, sys[2].position - sys[1].position);
            BOOST_TEST(dot1 == 0.0, boost::test_tools::tolerance(tol));
            BOOST_TEST(dot2 == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f1(0., 1., 0.);
            BOOST_TEST(mjolnir::math::length(normalize(sys[0].force) - f1) == 0.0, boost::test_tools::tolerance(tol));

            const coord_type f3(
                    cos(theta - mjolnir::math::constants<real_type>::pi / 2.),
                    sin(theta - mjolnir::math::constants<real_type>::pi / 2.),
                    0.);
            BOOST_TEST(mjolnir::math::length(normalize(sys[2].force) - f3) == 0.0, boost::test_tools::tolerance(tol));
        }

        // perpendicular to z axis
        BOOST_TEST(sys[0].force[2] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys[1].force[2] == 0.0, boost::test_tools::tolerance(tol));
        BOOST_TEST(sys[2].force[2] == 0.0, boost::test_tools::tolerance(tol));
    }
}

BOOST_AUTO_TEST_CASE(BondAngleInteraction_numerical_diff)
{
    using traits_type     = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type       = traits_type::real_type;
    using coord_type      = traits_type::coordinate_type;
    using boundary_type   = traits_type::boundary_type;
    using system_type     = mjolnir::System<traits_type>;
    using harmonic_type   = mjolnir::HarmonicPotential<real_type>;
    using bond_angle_type = mjolnir::BondAngleInteraction<traits_type, harmonic_type>;

    const real_type k(10.0);
    const real_type native(mjolnir::math::constants<real_type>::pi / 3.0);

    harmonic_type potential{k, native};
    bond_angle_type interaction("none", {{ {{0,1,2}}, potential}});

    system_type sys(3, boundary_type{});

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
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    constexpr real_type tol = 1e-5;
    constexpr real_type dr  = 1e-5;
    for(int i = 0; i < 1000; ++i)
    {
        for(std::size_t idx=0; idx<3; ++idx)
        {
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                const auto dx = uni(mt) * dr;

                mjolnir::math::X(sys.position(idx)) += dx;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::X(sys.position(idx)) += dx;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE == dx * mjolnir::math::X(sys.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                const auto dy = uni(mt) * dr;

                mjolnir::math::Y(sys.position(idx)) += dy;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::Y(sys.position(idx)) += dy;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE == dy * mjolnir::math::Y(sys.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                const auto dz = uni(mt) * dr;

                mjolnir::math::Z(sys.position(idx)) += dz;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::Z(sys.position(idx)) += dz;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE == dz * mjolnir::math::Z(sys.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
        }
    }
}


#include <mjolnir/interaction/local/bond_angle.hpp>

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
