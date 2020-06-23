#ifndef MJOLNIR_INPUT_READ_CONSTRAINT_INTERACTION_HPP
#define MJOLNIR_INPUT_READ_CONSTRAINT_INTERACTION_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/input/utility.hpp>
#include <mjolnir/core/ConstraintForceField.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

template<typename traitsT>
ConstraintForceField<traitsT>
read_constraint(const toml::value& constraint)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    using real_type            = typename traitsT::real_type;
    using indices_type         = std::array<std::size_t, 2>;
    using indices_v0_pair_type = std::pair<indices_type, real_type>;

    const auto kind          = toml::find<std::string>(constraint, "topology");
    const auto max_iteration = toml::find<std::size_t>(constraint, "max_iteration");
    const auto tolerance     = toml::find<real_type>  (constraint, "tolerance");
    const auto v0s  = toml::find<toml::array>(constraint, "v0s");
    MJOLNIR_LOG_NOTICE("-- ", v0s.size(), " constraints are found.");

    const auto& env = constraint.contains("env") ? constraint.at("env") : toml::value{};

    std::vector<indices_v0_pair_type> constraints;
    constraints.reserve(v0s.size());
    for(const auto& item : v0s)
    {
        const auto offset = find_parameter_or<std::int64_t>(item, env, "offset", 0);
        auto indices = find_parameter<indices_type>(item, env, "indices");
        auto v0      = find_parameter<real_type>   (item, env, "v0");
        for(auto& i : indices) {i += offset;}

        MJOLNIR_LOG_INFO_NO_LF("idxs = ", indices, ", v0 = ", v0);

        constraints.emplace_back(indices, v0);
    }
    // TODO : check parameter overlap in indices_v0s.
    return ConstraintForceField<traitsT>(kind, std::move(constraints),
            max_iteration, tolerance);
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template ConstraintForceField<SimulatorTraits<double, UnlimitedBoundary>>
  read_constraint(const toml::value& constraint);
extern template ConstraintForceField<SimulatorTraits<float,  UnlimitedBoundary>>
  read_constraint(const toml::value& constraint);
extern template ConstraintForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>>
  read_constraint(const toml::value& constraint);
extern template ConstraintForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>>
  read_constraint(const toml::value& constraint);
#endif

} // mjolnir
#endif /* MJOLNIR_INPUT_READ_CONSTRAINT_INTERACTION_HPP */
