#define BOOST_TEST_MODULE "test_afm_fitting_interaction"

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/AFMFit/AFMFitInteraction.hpp>
#include <mjolnir/util/make_unique.hpp>

#include <random>

namespace test
{
template<typename realT, typename coordT>
std::vector<realT>
make_image(const std::vector<realT>& radii, const std::vector<coordT>& positions,
           const realT pixel_width, const realT sigma, const realT gamma,
           const std::size_t n_x, const std::size_t n_y)
{
    std::vector<realT> image; image.reserve(n_x * n_y);
    for(std::size_t i_y=0; i_y < n_y; ++i_y)
    {
        const realT y0 = pixel_width * (i_y + 0.5);
        for(std::size_t i_x=0; i_x < n_x; ++i_x)
        {
            const realT x0 = pixel_width * (i_x + 0.5);

            realT sumexp = 1.0; // exp(0)
            for(std::size_t i=0; i<positions.size(); ++i)
            {
                const auto& p = positions.at(i);
                const auto& r = radii.at(i);
                const realT dx = mjolnir::math::X(p) - x0;
                const realT dy = mjolnir::math::Y(p) - y0;
                const realT dz = mjolnir::math::Z(p) - 0.0;

                const realT gx = -0.5 * dx * dx / (sigma * sigma);
                const realT gy = -0.5 * dy * dy / (sigma * sigma);
                const realT hz = (dz + r) / gamma;

                sumexp += std::exp(gx + gy + hz);
            }
            image.push_back(gamma * std::log(sumexp));
        }
    }
    assert(image.size() == n_x * n_y);
    return image;
}
}

BOOST_AUTO_TEST_CASE(AFMFitting_calc_force)
{
    mjolnir::LoggerManager::set_default_logger("test_afm_fitting_interaction.log");

    using traits_type      = mjolnir::SimulatorTraits<double, mjolnir::UnlimitedBoundary>;
    using real_type        = traits_type::real_type;
    using coord_type       = traits_type::coordinate_type;
    using boundary_type    = traits_type::boundary_type;
    using system_type      = mjolnir::System<traits_type>;
    using interaction_type = mjolnir::AFMFitInteraction<traits_type>;

    const std::size_t num_particles = 5;

    const real_type   k        =  5.0;
    const real_type   gamma    =  1.0;
    const real_type   z0       =  0.0;
    const real_type   cutoff   =  5.01;
    const real_type   margin   =  0.5;
    const real_type   sigma_x  =  2.0;
    const real_type   sigma_y  = sigma_x;
    const real_type   pixel_x  = 10.0;
    const real_type   pixel_y  = pixel_x;
    const std::size_t length_x = 10u;
    const std::size_t length_y = 10u;
    const std::vector<std::pair<std::size_t, real_type>> radii = {
        {1, 1.0}, {2, 2.0}, {3, 3.0}, {4, 4.0}, {5, 5.0}
    };

    std::vector<coord_type> initial_coordinates = {
        coord_type(20.0, 20.0, 10.0),
        coord_type(50.0, 50.0, 15.0),
        coord_type(80.0, 80.0,  5.0),
        coord_type(20.0, 50.0, 20.0),
        coord_type(50.0, 80.0, 10.0)
    };

    const std::vector<real_type> image = test::make_image(
            {1.0, 2.0, 3.0, 4.0, 5.0}, initial_coordinates,
            pixel_x, sigma_x, gamma, length_x, length_y);

    auto& l_o_g_g_e_r_ = mjolnir::LoggerManager::get_default_logger();
    l_o_g_g_e_r_.log(mjolnir::Logger::Level::Notice, image);

    interaction_type interaction(k, gamma, z0, cutoff, margin,
            sigma_x, sigma_y, pixel_x, pixel_y, length_x, length_y,
            radii, image);

    std::mt19937 mt(123456789);
    std::uniform_real_distribution<real_type> uni(-1.0, 1.0);

    for(std::size_t trial = 0; trial < 1000; ++trial)
    {
        system_type sys(num_particles, boundary_type{});

        for(std::size_t i=0; i<sys.size(); ++i)
        {
            sys.mass(i)  = 1.0;
            sys.rmass(i) = 1.0;

            sys.position(i)  = initial_coordinates.at(i);
            sys.position(i) += coord_type(0.01 * uni(mt), 0.01 * uni(mt), 0.01 * uni(mt));
            sys.velocity(i)  = coord_type(0.0, 0.0, 0.0);
            sys.force(i)     = coord_type(0.0, 0.0, 0.0);
            sys.name(i)      = "X";
            sys.group(i)     = "TEST";
        }

        const auto init = sys;
        interaction.initialize(init);

        constexpr real_type tol = 1e-3;
        constexpr real_type dr  = 1e-4;
        for(std::size_t idx=0; idx<sys.size(); ++idx)
        {
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                mjolnir::math::X(sys.position(idx)) += dr;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::X(sys.position(idx)) += dr;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE == dr * mjolnir::math::X(sys.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                mjolnir::math::Y(sys.position(idx)) += dr;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::Y(sys.position(idx)) += dr;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE == dr * mjolnir::math::Y(sys.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
            {
                // ----------------------------------------------------------------
                // reset positions
                sys = init;

                // calc U(x-dx)
                const auto E0 = interaction.calc_energy(sys);

                mjolnir::math::Z(sys.position(idx)) += dr;

                // calc F(x)
                interaction.calc_force(sys);

                mjolnir::math::Z(sys.position(idx)) += dr;

                // calc U(x+dx)
                const auto E1 = interaction.calc_energy(sys);

                // central difference
                const auto dE = (E1 - E0) * 0.5;

                BOOST_TEST(-dE == dr * mjolnir::math::Z(sys.force(idx)),
                           boost::test_tools::tolerance(tol));
            }
        }
    }
}
