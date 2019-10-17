#ifndef MJOLNIR_INPUT_READ_LOCAL_INTERACTION_HPP
#define MJOLNIR_INPUT_READ_LOCAL_INTERACTION_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/interaction/local/BondLengthInteraction.hpp>
#include <mjolnir/interaction/local/ContactInteraction.hpp>
#include <mjolnir/interaction/local/DirectionalContactInteraction.hpp>
#include <mjolnir/interaction/local/BondAngleInteraction.hpp>
#include <mjolnir/interaction/local/DihedralAngleInteraction.hpp>
#include <mjolnir/interaction/local/DummyInteraction.hpp>

#include <mjolnir/forcefield/3SPN2/ThreeSPN2BaseStackingPotential.hpp>
#include <mjolnir/forcefield/3SPN2/ThreeSPN2BaseStackingInteraction.hpp>
#include <mjolnir/forcefield/PDNS/ProteinDNANonSpecificInteraction.hpp>

#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_local_potential.hpp>
#include <memory>
#include <utility>

namespace mjolnir
{

// ----------------------------------------------------------------------------
// individual local interaction. would be called from read_local_interaction
// defined at the bottom of this file.
// ----------------------------------------------------------------------------

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_bond_length_interaction(const std::string& kind, const toml::value& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = typename traitsT::real_type;

    const auto potential = toml::find<std::string>(local, "potential");
    if(potential == "Harmonic")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Harmonic.");
        using potentialT = HarmonicPotential<real_type>;

        return make_unique<BondLengthInteraction<traitsT, potentialT>>(
                kind, read_local_potential<2, potentialT>(local));
    }
    else if(potential == "GoContact")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is 10-12 Go contact.");
        using potentialT = GoContactPotential<real_type>;

        return make_unique<BondLengthInteraction<traitsT, potentialT>>(
                kind, read_local_potential<2, potentialT>(local));
    }
    else if(potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Gaussian.");
        using potentialT = GaussianPotential<real_type>;

        return make_unique<BondLengthInteraction<traitsT, potentialT>>(
                kind, read_local_potential<2, potentialT>(local));
    }
    else if(potential == "3SPN2Bond")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is 3SPN2Bond.");
        using potentialT = ThreeSPN2BondPotential<real_type>;

        return make_unique<BondLengthInteraction<traitsT, potentialT>>(
                kind, read_local_potential<2, potentialT>(local));
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_bond_length_interaction: invalid potential",
            toml::find<toml::value>(local, "potential"), "here", {
            "expected value is one of the following.",
            "- \"Harmonic\" : well-known harmonic potential",
            "- \"GoContact\": r^12 - r^10 type native contact potential",
            "- \"Gaussian\" : well-known gaussian potential",
            "- \"3SPN2Bond\": bond length potential for 3SPN2"
            }));
    }
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_contact_interaction(const std::string& kind, const toml::value& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = typename traitsT::real_type;

    real_type margin = 0.5; // default value
    if(local.as_table().count("margin") == 1)
    {
        margin = toml::find<real_type>(local, "margin");
    }

    const auto potential = toml::find<std::string>(local, "potential");
    if(potential == "GoContact")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is 10-12 Go contact.");
        using potentialT = GoContactPotential<real_type>;

        return make_unique<ContactInteraction<traitsT, potentialT>>(
                kind, read_local_potential<2, potentialT>(local), margin);
    }
    else if(potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Gaussian.");
        using potentialT = GaussianPotential<real_type>;

        return make_unique<ContactInteraction<traitsT, potentialT>>(
                kind, read_local_potential<2, potentialT>(local), margin);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_bond_length_interaction: invalid potential",
            toml::find<toml::value>(local, "potential"), "here", {
            "expected value is one of the following.",
            "- \"GoContact\": r^12 - r^10 type native contact potential",
            "- \"Gaussian\" : well-known gaussian potential"
            }));
    }
}

template<std::size_t N, typename realT,
         typename angle1T, typename angle2T, typename contactT>
std::vector<std::tuple<std::array<std::size_t, N>, angle1T, angle2T, contactT>>
read_directional_contact_potentials(const toml::value& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    MJOLNIR_LOG_INFO("as", N, "-body interaction");

    using indices_type = std::array<std::size_t, N>;
    using indices_potentials_tuple_type =
        std::tuple<indices_type, angle1T, angle2T, contactT>;

    const auto& params = toml::find<toml::array>(local, "parameters");
    MJOLNIR_LOG_NOTICE("-- ", params.size(), " interactions are found.");

    const auto& env = local.as_table().count("env") == 1 ?
                      local.as_table().at("env") : toml::value{};

    std::vector<indices_potentials_tuple_type> retval;
    retval.reserve(params.size());
    for(const auto& item : params)
    {
        const auto indices = find_parameter<indices_type>(item, env, "indices");
        MJOLNIR_LOG_INFO_NO_LF("idxs = ", indices, ", ");

        const auto angle1  = find_parameter<toml::value>(item, env, "angle1");
        const auto angle2  = find_parameter<toml::value>(item, env, "angle2");
        const auto contact = find_parameter<toml::value>(item, env, "contact");

        retval.push_back(std::make_tuple(indices,
            detail::read_local_potential_impl<angle1T>::invoke(angle1, env),
            detail::read_local_potential_impl<angle2T>::invoke(angle2, env),
            detail::read_local_potential_impl<contactT>::invoke(contact, env)));
    }
    return retval;
}

// This is reading contact part of read_directional_contact_interaction function.
template<typename traitsT, typename ... PotentialTs>
typename std::enable_if<sizeof...(PotentialTs) == 2,
    std::unique_ptr<LocalInteractionBase<traitsT>>>::type
read_directional_contact_interaction(const std::string& kind,
        const toml::value& local, std::vector<std::string>)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = typename traitsT::real_type;

    real_type mgn = 0.5; // default value
    if(local.as_table().count("margin") == 1)
    {
        mgn = toml::find<real_type>(local, "margin");
    }
    const real_type margin = mgn;

    const auto contact_potential = toml::find<std::string>(local, "potentials", "contact");
    if(contact_potential == "GoContact")
    {
        MJOLNIR_LOG_NOTICE("-- contact potential function is 10-12 Go contact.");
        using contact_potentialT = GoContactPotential<real_type>;

        return make_unique<DirectionalContactInteraction<
            traitsT, PotentialTs..., contact_potentialT>>(kind,
              read_directional_contact_potentials<
                  4, real_type, PotentialTs..., contact_potentialT
              >(local), margin);
    }
    else if(contact_potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- contact potential function is Gaussian.");
        using contact_potentialT = GaussianPotential<real_type>;

        return make_unique<DirectionalContactInteraction<
                traitsT, PotentialTs..., contact_potentialT>>(kind,
              read_directional_contact_potentials<
                  4, real_type, PotentialTs..., contact_potentialT
              >(local), margin);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_directional_contact_interaction: invalid contact potential",
            toml::find<toml::value>(local, "potentials", "contact"), "here", {
            "expected value is one of the following.",
            "- \"GoContact\": r^12 - r^10 type native contact potential",
            "- \"Gaussian\" : well-known gaussian potential"
            }));
    }
}

// This is reading angle parts of read_directional_contact_interaction function.
template<typename traitsT, typename ... PotentialTs>
typename std::enable_if<sizeof...(PotentialTs) < 2,
    std::unique_ptr<LocalInteractionBase<traitsT>>>::type
read_directional_contact_interaction(const std::string& kind,
        const toml::value& local, std::vector<std::string> angle_keys)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    using real_type = typename traitsT::real_type;

    const auto angle_key = angle_keys.back();
    const auto angle_potential =
        toml::find<std::string>(local, "potentials", angle_key);

    angle_keys.pop_back();
    if(angle_potential == "Cosine")
    {
        MJOLNIR_LOG_NOTICE("-- angle potential function is Cosine potential");
        using angle_potential_T = CosinePotential<real_type>;

        return read_directional_contact_interaction<
            traitsT, PotentialTs..., angle_potential_T
            >(kind, local, std::move(angle_keys));
    }
    else if(angle_potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- angle potential function is Gaussian");
        using angle_potential_T = GaussianPotential<real_type>;

        return read_directional_contact_interaction<
            traitsT, PotentialTs..., angle_potential_T
            >(kind, local, std::move(angle_keys));
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_bond_length_interaction: invalid angle potential",
            toml::find(local, "potentials", angle_key), "here", {
            "expected value is one of the following.",
            "- \"Cosine\"   : 1 + Cosine(x) potential",
            "- \"Gaussian\" : well-known gaussian potential"
            }));
    }
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_bond_angle_interaction(const std::string& kind, const toml::value& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = typename traitsT::real_type;

    const auto potential = toml::find<std::string>(local, "potential");

    if(potential == "Harmonic")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Harmonic.");
        using potentialT = HarmonicPotential<real_type>;

        return make_unique<BondAngleInteraction<traitsT, potentialT>>(
                kind, read_local_potential<3, potentialT>(local));
    }
    else if(potential == "FlexibleLocalAngle")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Flexible Local Angle.");
        using potentialT = FlexibleLocalAnglePotential<real_type>;

        return make_unique<BondAngleInteraction<traitsT, potentialT>>(
                kind, read_local_potential<3, potentialT>(local));
    }
    else if(potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Gaussian.");
        using potentialT = GaussianPotential<real_type>;

        return make_unique<BondAngleInteraction<traitsT, potentialT>>(
                kind, read_local_potential<3, potentialT>(local));
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_bond_angle_interaction: invalid potential",
            toml::find<toml::value>(local, "potential"), "here", {
            "expected value is one of the following.",
            "- \"Harmonic\"          : well-known harmonic potential",
            "- \"Gaussian\"          : well-known gaussian potential"
            "- \"FlexibleLocalAngle\": table-based potential for C-alpha protein model",
            }));
    }
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_dihedral_angle_interaction(const std::string& kind, const toml::value& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type = typename traitsT::real_type;

    const auto potential = toml::find<std::string>(local, "potential");

    if(potential == "ClementiDihedral")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Clementi-Go's dihedral.");
        using potentialT = ClementiDihedralPotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    else if(potential == "Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Gaussian.");
        using potentialT = PeriodicGaussianPotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    else if(potential == "FlexibleLocalDihedral")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Flexible Local Dihedral.");
        using potentialT = FlexibleLocalDihedralPotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    else if(potential == "Cosine")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is Cosine.");
        using potentialT = CosinePotential<real_type>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(
            kind, read_local_potential<4, potentialT>(local));
    }
    else if(potential == "Gaussian+FlexibleLocalDihedral" ||
            potential == "FlexibleLocalDihedral+Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is "
                           "Gaussian + FlexibleLocalDihedral.");
        using potentialT = SumLocalPotential<real_type,
              PeriodicGaussianPotential, FlexibleLocalDihedralPotential>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(kind,
            read_local_potentials<4, real_type,
                PeriodicGaussianPotential, FlexibleLocalDihedralPotential
                >(local, "Gaussian", "FlexibleLocalDihedral"));
    }
    else if(potential == "Gaussian+Cosine" || potential == "Cosine+Gaussian")
    {
        MJOLNIR_LOG_NOTICE("-- potential function is PeriodicGaussian + Cosine.");
        using potentialT = SumLocalPotential<real_type,
              PeriodicGaussianPotential, CosinePotential>;

        return make_unique<DihedralAngleInteraction<traitsT, potentialT>>(kind,
            read_local_potentials<
                4, real_type, PeriodicGaussianPotential, CosinePotential
                >(local, "Gaussian", "Cosine"));
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_dihedral_angle_interaction: invalid potential",
            toml::find<toml::value>(local, "potential"), "here", {
            "expected value is one of the following.",
            "- \"Gaussian\"             : well-known gaussian potential"
            "- \"ClementiDihedral\"     : potential used in the off-lattice Go protein model"
            "- \"FlexibleLocalDihedral\": table-based potential for C-alpha protein model",
            }));
    }
}


template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_dummy_interaction(const std::string& kind, const toml::value& local)
{
    // It does not require `potential` field because this interaction does not
    // calculate force or energy.
    //
    // ```toml
    // [[forcefields.local]]
    // interaction = "Dummy"
    // topology  = "bond"
    // parameters = [
    //     {indices = [0, 1]},
    //     # ...
    // ]
    // ```
    //
    // So, unlike other interactions, it reads particle indices here.

    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using indices_type = std::array<std::size_t, 2>;

    if(kind == "none")
    {
        MJOLNIR_LOG_WARN("Dummy potential got `none`. It has no effect.");
        MJOLNIR_LOG_WARN("It is for defining unusual topology.");
    }
    if(local.as_table().count("potential") == 1)
    {
        MJOLNIR_LOG_WARN("Dummy Interaction has a `potential` field.");
        MJOLNIR_LOG_WARN("It is for defining unusual topology.");
        MJOLNIR_LOG_WARN("Potential parameters are ignored.");
    }

    const auto& params = toml::find<toml::array>(local, "parameters");
    MJOLNIR_LOG_NOTICE("-- ", params.size(), " bonds are found.");

    std::vector<indices_type> indices_list;
    indices_list.reserve(params.size());
    for(const auto& item : params)
    {
        const auto indices = toml::find<indices_type>(item, "indices");
        MJOLNIR_LOG_INFO("idxs = ", indices);
        indices_list.push_back(indices);
    }
    return make_unique<DummyInteraction<traitsT>>(kind, std::move(indices_list));
}

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_3spn2_base_stacking_interaction(const std::string& kind, const toml::value& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type      = typename traitsT::real_type;
    using potential_type = ThreeSPN2BaseStackingPotential<real_type>;
    using base_kind      = parameter_3SPN2::base_kind;

    // [[forcefields.local]]
    // interaction = "3SPN2BaseStacking"
    // topology    = "nucleotide"
    // parameters = [
    //     {strand = 0, nucleotide =  0,          S =   0, B =   1, Base = "A"},
    //     {strand = 0, nucleotide =  1, P =   2, S =   3, B =   4, Base = "T"},
    //     # ...
    // ]

    using indices_type   = std::array<std::size_t, 3>; // {S, B3, B5}
    using parameter_type = typename potential_type::parameter_type;

    potential_type potential;
    const auto pot = toml::find<std::string>(local, "potential");
    if(pot == "3SPN2")
    {
        MJOLNIR_LOG_NOTICE("potential is 3SPN2");
        ThreeSPN2BaseStackingPotentialParameter<real_type> para_3SPN2;
        potential = potential_type(para_3SPN2);
    }
    else if (pot == "3SPN2C")
    {
        MJOLNIR_LOG_NOTICE("potential is 3SPN2C");
        ThreeSPN2CBaseStackingPotentialParameter<real_type> para_3SPN2C;
        potential = potential_type(para_3SPN2C);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_local_3spn2_base_stacking_interaction: "
            "invalid potential", toml::find(local, "potential"), "here", {
            "expected value is one of the following.",
            "- \"3SPN2\" : The general 3SPN2 parameter set.",
            "- \"3SPN2C\": The parameter set optimized to reproduce curveture of dsDNA."
            }));
    }

    const auto& params = toml::find<toml::array>(local, "parameters");
    MJOLNIR_LOG_NOTICE("-- ", params.size(), " interactions are found.");

    const auto& env = local.as_table().count("env") == 1 ?
                      local.as_table().at("env") : toml::value{};
    // parameters = [
    //     {strand = 0, nucleotide =  0,          S =   0, B =   1, Base = "A"},
    //     {strand = 0, nucleotide =  1, P =   2, S =   3, B =   4, Base = "T"},
    //     # ...
    // ]
    using nucleotide_index_type = parameter_3SPN2::NucleotideIndex;
    std::vector<nucleotide_index_type> nuc_idxs;
    nuc_idxs.reserve(params.size());
    for(const auto& item : params)
    {
        nucleotide_index_type nuc_idx;

        // at the edge of the DNA, Phosphate may not exist.
        if(item.as_table().count("P") != 0)
        {
            nuc_idx.P = find_parameter<std::size_t>(item, env, "P");
        }
        nuc_idx.S          = find_parameter<std::size_t>(item, env, "S");
        nuc_idx.B          = find_parameter<std::size_t>(item, env, "B");
        nuc_idx.strand     = find_parameter<std::size_t>(item, env, "strand");
        nuc_idx.nucleotide = find_parameter<std::size_t>(item, env, "nucleotide");

        const auto bk      = find_parameter<std::string>(item, env, "Base");
        if     (bk == "A") {nuc_idx.base = base_kind::A;}
        else if(bk == "T") {nuc_idx.base = base_kind::T;}
        else if(bk == "G") {nuc_idx.base = base_kind::G;}
        else if(bk == "C") {nuc_idx.base = base_kind::C;}
        else
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_local_3spn2_base_stacking_interaction: "
                "invalid Base", item, "here", {
                "expected value is one of the \"A\", \"T\", \"C\", \"G\"."
                }));
        }

        MJOLNIR_LOG_INFO("ThreeSPN2BaseStackingPotential: nucleotide = ", nuc_idx);
        nuc_idxs.push_back(nuc_idx);
    }

    std::sort(nuc_idxs.begin(), nuc_idxs.end(),
        [](const nucleotide_index_type& lhs, const nucleotide_index_type& rhs) {
            return std::make_pair(lhs.strand, lhs.nucleotide) <
                   std::make_pair(rhs.strand, rhs.nucleotide);
        });

    std::vector<std::pair<indices_type, parameter_type>> parameters;
    parameters.reserve(params.size());
    for(std::size_t i=1; i<nuc_idxs.size(); ++i)
    {
        const auto& Base5 = nuc_idxs.at(i-1);
        const auto& Base3 = nuc_idxs.at(i);

        if(Base5.strand != Base3.strand)
        {
            continue; // if the strands are different, there is no stacking.
        }
        assert(Base3.base != base_kind::X);
        assert(Base5.base != base_kind::X);

        const std::array<std::size_t, 3> idxs{{Base5.S, Base5.B, Base3.B}};
        const auto bs_kind = potential.bs_kind(Base5.base, Base3.base);

        MJOLNIR_LOG_INFO("ThreeSPN2BaseStackingPotential = {S = ", Base5.S,
            ", B5 = ", Base5.B, ", B3 = ", Base3.B, ", bases = ", bs_kind, '}');

        parameters.emplace_back(idxs, bs_kind);
    }
    return make_unique<ThreeSPN2BaseStackingInteraction<traitsT>>(kind,
            std::move(parameters), std::move(potential), std::move(nuc_idxs));
}


// ----------------------------------------------------------------------------
// PDNS Interaction

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_pdns_interaction(const std::string& topol, const toml::value& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type      = typename traitsT::real_type;
    using potential_type = ProteinDNANonSpecificPotential<real_type>;
    using parameter_type = typename potential_type::parameter_type;
    using dna_index_type = typename potential_type::dna_index_type;

    if(topol != "none")
    {
        MJOLNIR_LOG_WARN("In PDNS Interaction, contacting particle pair changes"
                         " every step. Topology cannot be defined.");
    }
    // ```toml
    // [[forcefields.local]]
    // interaction = "PDNS"
    // potential   = "PDNS"
    // topology    = "none"
    // sigma  = 1.0
    // delta  = 0.17453
    // cutoff = 5.0 # relative to sigma
    // parameters  = [
    // {index =    2, kind = "DNA", S3 = 1},
    // {index =    5, kind = "DNA", S3 = 4},
    // # ...
    // {index = 1000, kind = "Protein", PN =  999, PC = 1001, k = 1.2, r0 = 5.0, theta0 = 100.0, phi0 = 130.0},
    // {index = 1023, kind = "Protein", PN = 1022, PC = 1024, k = 1.2, r0 = 6.0, theta0 = 110.0, phi0 = 120.0},
    // # ...
    // ]
    // ```

    // ------------------------------------------------------------------------
    // read parameters

    const real_type dlt = toml::find<real_type>(local, "delta");
    const real_type sgm = toml::find<real_type>(local, "sigma");

    real_type cut = potential_type::default_cutoff();
    if(local.as_table().count("cutoff") != 0)
    {
        cut = toml::find<real_type>(local, "cutoff");
    }
    const real_type cutoff = cut;
    const real_type margin = toml::find_or<real_type>(local, "margin", 0.5);

    const auto& env = local.as_table().count("env") == 1 ?
                      local.as_table().at("env") : toml::value{};

    const auto& ps = toml::find<toml::array>(local, "parameters");
    MJOLNIR_LOG_INFO(ps.size(), " parameters are found");

    std::vector<parameter_type> params;
    std::vector<dna_index_type> dnas;
    params.reserve(ps.size());
    dnas  .reserve(ps.size());

    for(const auto& item : ps)
    {
        const auto idx  = toml::find<std::size_t>(item, "index");
        const auto kind = toml::find<std::string>(item, "kind");
        if(kind == "Protein")
        {
            parameter_type para;
            para.P      = idx;
            para.PN     = find_parameter<std::uint32_t>(item, env, "PN");
            para.PC     = find_parameter<std::uint32_t>(item, env, "PC");
            para.k      = find_parameter<real_type>(item, env, "k");
            para.r0     = find_parameter<real_type>(item, env, "r0");
            para.theta0 = find_parameter<real_type>(item, env, "theta0");
            para.phi0   = find_parameter<real_type>(item, env, "phi0");
            params.push_back(para);

            MJOLNIR_LOG_INFO("Protein: idx = ", idx, ", PN = ", para.PN,
                ", PC = ", para.PC, ", k = ", para.k, ", r0 = ", para.r0,
                ", theta0 = ", para.theta0, ", phi0 = ", para.phi0);
        }
        else if (kind == "DNA")
        {
            const auto S3 = find_parameter<std::uint32_t>(item, env, "S3");
            dnas.emplace_back(idx, S3);
            MJOLNIR_LOG_INFO("DNA: idx = ", idx, ", S3 = ", S3);
        }
        else
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_pdns_interaction: unknown kind ",
                toml::find(item, "kind"), "here", {
                "expected value is one of the following.",
                "- \"Protein\": Protein bead has PN, PC and native parameters",
                "- \"DNA\"    : DNA bead has index of the corresponding Sugar."
                }));
        }
    }
    return make_unique<ProteinDNANonSpecificInteraction<traitsT>>(
        potential_type(sgm, dlt, cutoff, std::move(params), std::move(dnas)),
        margin);
}

// ----------------------------------------------------------------------------
// general read_local_interaction function
// ----------------------------------------------------------------------------

template<typename traitsT>
std::unique_ptr<LocalInteractionBase<traitsT>>
read_local_interaction(const toml::value& local)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto interaction = toml::find<std::string>(local, "interaction");
    const auto kind        = toml::find<std::string>(local, "topology");
    MJOLNIR_LOG_INFO("connection kind = ", kind);

    if(interaction == "BondLength")
    {
        MJOLNIR_LOG_NOTICE("Bond Length interaction found.");
        return read_bond_length_interaction<traitsT>(kind, local);
    }
    else if(interaction == "Contact")
    {
        MJOLNIR_LOG_NOTICE("Contact interaction found.");
        return read_contact_interaction<traitsT>(kind, local);
    }
    else if(interaction == "DirectionalContact")
    {
        MJOLNIR_LOG_NOTICE("Directional Contact interaction found.");
        return read_directional_contact_interaction<traitsT>(
                kind, local, {"angle2", "angle1"});
    }
    else if(interaction == "BondAngle")
    {
        MJOLNIR_LOG_NOTICE("Bond Angle interaction found.");
        return read_bond_angle_interaction<traitsT>(kind, local);
    }
    else if(interaction == "DihedralAngle")
    {
        MJOLNIR_LOG_NOTICE("Dihedral Angle interaction found.");
        return read_dihedral_angle_interaction<traitsT>(kind, local);
    }
    else if(interaction == "Dummy")
    {
        MJOLNIR_LOG_NOTICE("Dummy interaction found.");
        return read_dummy_interaction<traitsT>(kind, local);
    }
    else if(interaction == "3SPN2BaseStacking")
    {
        MJOLNIR_LOG_NOTICE("3SPN2BaseStacking interaction found.");
        return read_3spn2_base_stacking_interaction<traitsT>(kind, local);
    }
    else if(interaction == "PDNS")
    {
        MJOLNIR_LOG_NOTICE("P-D ns Interaction found.");
        return read_pdns_interaction<traitsT>(kind, local);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_local_interaction: invalid interaction",
            toml::find<toml::value>(local, "interaction"), "here", {
            "expected value is one of the following.",
            "- \"BondLength\"    : 2-body well-known chemical bond interaction",
            "- \"BondAngle\"     : 3-body well-known bond angle interaction",
            "- \"Contact\"       : 2-body bond interaction that might be broken",
            "- \"DihedralAngle\" : 4-body well-known dihedral angle interaction",
            "- \"DirectionalContact\" : 4-body contact interaction depends on the contact angle",
            "- \"Dummy\"         : To represent a strange topology. It does nothing",
            "- \"PDNS\"          : directional contact representing H-bond between protein and DNA"
            }));
    }
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_local_interaction(const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_local_interaction(const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_local_interaction(const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_local_interaction(const toml::value& local);

extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_bond_length_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_bond_length_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_bond_length_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_bond_length_interaction(const std::string& kind, const toml::value& local);

extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_contact_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_contact_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_contact_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_contact_interaction(const std::string& kind, const toml::value& local);

extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_bond_angle_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_bond_angle_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_bond_angle_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_bond_angle_interaction(const std::string& kind, const toml::value& local);

extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_dihedral_angle_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_dihedral_angle_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_dihedral_angle_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_dihedral_angle_interaction(const std::string& kind, const toml::value& local);

extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_dummy_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_dummy_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_dummy_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_dummy_interaction(const std::string& kind, const toml::value& local);

extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, UnlimitedBoundary>       >> read_3spn2_base_stacking_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  UnlimitedBoundary>       >> read_3spn2_base_stacking_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>> read_3spn2_base_stacking_interaction(const std::string& kind, const toml::value& local);
extern template std::unique_ptr<LocalInteractionBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>> read_3spn2_base_stacking_interaction(const std::string& kind, const toml::value& local);
#endif


} // mjolnir
#endif// MJOLNIR_READ_LOCAL_INTERACTION
