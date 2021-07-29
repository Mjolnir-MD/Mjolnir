#ifndef MJOLNIR_TEST_OMP_CHECK_CONSISTENCY_HPP
#define MJOLNIR_TEST_OMP_CHECK_CONSISTENCY_HPP
#include <mjolnir/math/Vector.hpp>
#include <mjolnir/core/System.hpp>

#include <test/util/clear_system.hpp>

#ifdef BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#else
#include <boost/test/included/unit_test.hpp>
#endif

namespace mjolnir
{
namespace test
{

template<typename SequentialTrait, typename SequentialInteraction,
         typename ParallelTrait,   typename ParallelInteraction>
void check_force_consistency(System<SequentialTrait>      seq_sys,
                             const SequentialInteraction& seq_interaction,
                             System<ParallelTrait>        par_sys,
                             const ParallelInteraction    par_interaction,
                             const typename SequentialTrait::real_type tol)
{
    BOOST_TEST_REQUIRE(seq_sys.size() == par_sys.size());

    clear_force(seq_sys);
    seq_sys.preprocess_forces();
    seq_interaction.calc_force(seq_sys);
    seq_sys.postprocess_forces();

    clear_force(par_sys);
    par_sys.preprocess_forces();
    par_interaction.calc_force(par_sys);
    par_sys.postprocess_forces();

    // check the values are the same
    for(std::size_t i=0; i<seq_sys.size(); ++i)
    {
        BOOST_TEST(mjolnir::math::X(seq_sys.force(i)) == mjolnir::math::X(par_sys.force(i)),
                   boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Y(seq_sys.force(i)) == mjolnir::math::Y(par_sys.force(i)),
                   boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Z(seq_sys.force(i)) == mjolnir::math::Z(par_sys.force(i)),
                   boost::test_tools::tolerance(tol));
    }
    // check the virials are the same
    for(std::size_t i=0; i<9; ++i)
    {
        BOOST_TEST(seq_sys.virial()[i] == par_sys.virial()[i], boost::test_tools::tolerance(tol));
    }

    return;
}

template<typename SequentialTrait, typename SequentialInteraction,
         typename ParallelTrait,   typename ParallelInteraction>
void check_force_and_virial_consistency(
        System<SequentialTrait>      seq_sys,
        const SequentialInteraction& seq_interaction,
        System<ParallelTrait>        par_sys,
        const ParallelInteraction    par_interaction,
        const typename SequentialTrait::real_type tol)
{
    BOOST_TEST_REQUIRE(seq_sys.size() == par_sys.size());

    clear_force(seq_sys);
    seq_sys.preprocess_forces();
    seq_interaction.calc_force_and_virial(seq_sys);
    seq_sys.postprocess_forces();

    clear_force(par_sys);
    par_sys.preprocess_forces();
    par_interaction.calc_force_and_virial(par_sys);
    par_sys.postprocess_forces();

    // check the values are the same
    for(std::size_t i=0; i<seq_sys.size(); ++i)
    {
        BOOST_TEST(mjolnir::math::X(seq_sys.force(i)) == mjolnir::math::X(par_sys.force(i)),
                   boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Y(seq_sys.force(i)) == mjolnir::math::Y(par_sys.force(i)),
                   boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Z(seq_sys.force(i)) == mjolnir::math::Z(par_sys.force(i)),
                   boost::test_tools::tolerance(tol));
    }
    // check the virials are the same
    for(std::size_t i=0; i<9; ++i)
    {
        BOOST_TEST(seq_sys.virial()[i] == par_sys.virial()[i], boost::test_tools::tolerance(tol));
    }

    return;
}


template<typename SequentialTrait, typename SequentialInteraction,
         typename ParallelTrait,   typename ParallelInteraction>
void check_force_and_energy_consistency(
        System<SequentialTrait>      seq_sys,
        const SequentialInteraction& seq_interaction,
        System<ParallelTrait>        par_sys,
        const ParallelInteraction    par_interaction,
        const typename SequentialTrait::real_type tol)
{
    BOOST_TEST_REQUIRE(seq_sys.size() == par_sys.size());

    clear_force(seq_sys);
    seq_sys.preprocess_forces();
    const auto seq_ene = seq_interaction.calc_force_and_energy(seq_sys);
    seq_sys.postprocess_forces();

    clear_force(par_sys);
    par_sys.preprocess_forces();
    const auto par_ene = par_interaction.calc_force_and_energy(par_sys);
    par_sys.postprocess_forces();

    BOOST_TEST(seq_ene == par_ene, boost::test_tools::tolerance(tol));

    // check the values are the same
    for(std::size_t i=0; i<seq_sys.size(); ++i)
    {
        BOOST_TEST(mjolnir::math::X(seq_sys.force(i)) == mjolnir::math::X(par_sys.force(i)),
                   boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Y(seq_sys.force(i)) == mjolnir::math::Y(par_sys.force(i)),
                   boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Z(seq_sys.force(i)) == mjolnir::math::Z(par_sys.force(i)),
                   boost::test_tools::tolerance(tol));
    }
    // check the virials are the same
    for(std::size_t i=0; i<9; ++i)
    {
        BOOST_TEST(seq_sys.virial()[i] == par_sys.virial()[i], boost::test_tools::tolerance(tol));
    }

    return;
}

template<typename SequentialTrait, typename SequentialInteraction,
         typename ParallelTrait,   typename ParallelInteraction>
void check_force_energy_virial_consistency(
        System<SequentialTrait>      seq_sys,
        const SequentialInteraction& seq_interaction,
        System<ParallelTrait>        par_sys,
        const ParallelInteraction    par_interaction,
        const typename SequentialTrait::real_type tol)
{
    clear_force(seq_sys);
    seq_sys.preprocess_forces();
    const auto seq_ene = seq_interaction.calc_force_virial_energy(seq_sys);
    seq_sys.postprocess_forces();

    clear_force(par_sys);
    par_sys.preprocess_forces();
    const auto par_ene = par_interaction.calc_force_virial_energy(par_sys);
    par_sys.postprocess_forces();

    BOOST_TEST(seq_ene == par_ene, boost::test_tools::tolerance(tol));

    // check the values are the same
    for(std::size_t i=0; i<seq_sys.size(); ++i)
    {
        BOOST_TEST(mjolnir::math::X(seq_sys.force(i)) == mjolnir::math::X(par_sys.force(i)),
                   boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Y(seq_sys.force(i)) == mjolnir::math::Y(par_sys.force(i)),
                   boost::test_tools::tolerance(tol));
        BOOST_TEST(mjolnir::math::Z(seq_sys.force(i)) == mjolnir::math::Z(par_sys.force(i)),
                   boost::test_tools::tolerance(tol));
    }
    // check the virials are the same
    for(std::size_t i=0; i<9; ++i)
    {
        BOOST_TEST(seq_sys.virial()[i] == par_sys.virial()[i], boost::test_tools::tolerance(tol));
    }

    return;
}

template<typename SequentialTrait, typename SequentialInteraction,
         typename ParallelTrait,   typename ParallelInteraction>
void check_energy_consistency(
        const System<SequentialTrait>& seq_sys,
        const SequentialInteraction&   seq_interaction,
        const System<ParallelTrait>&   par_sys,
        const ParallelInteraction      par_interaction,
        const typename SequentialTrait::real_type tol)
{
    const auto seq_ene = seq_interaction.calc_energy(seq_sys);
    const auto par_ene = par_interaction.calc_energy(par_sys);

    BOOST_TEST(seq_ene == par_ene, boost::test_tools::tolerance(tol));

    return;
}

} // test
} // mjolnir
#endif
