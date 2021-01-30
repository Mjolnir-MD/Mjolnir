#ifndef MJOLNIR_CORE_ZORDER_R_TREE_HPP
#define MJOLNIR_CORE_ZORDER_R_TREE_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/SpatialPartitionBase.hpp>
#include <mjolnir/util/range.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/fixed_vector.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

namespace mjolnir
{

template<typename traitsT, typename PotentialT, std::size_t MaxElem = 8>
class ZorderRTree final : public SpatialPartitionBase<traitsT, PotentialT>
{
    static_assert(MaxElem > 1, "MaxElem > 1");

  public:
    using traits_type        = traitsT;
    using potential_type     = PotentialT;
    using base_type          = SpatialPartitionBase<traits_type, potential_type>;

    using system_type        = typename base_type::system_type;
    using boundary_type      = typename base_type::boundary_type;
    using real_type          = typename base_type::real_type;
    using coordinate_type    = typename base_type::coordinate_type;
    using neighbor_list_type = typename base_type::neighbor_list_type;
    using neighbor_type      = typename base_type::neighbor_type;
    using range_type         = typename base_type::range_type;

    struct AABB
    {
        coordinate_type lower;
        coordinate_type upper;
    };
    struct Node
    {
        fixed_vector<std::size_t, MaxElem> children;
        std::size_t parent;
        AABB box;
        bool is_leaf;
    };

  public:

    ZorderRTree(): cutoff_(0), margin_(1), current_margin_(-1) {}
    ~ZorderRTree() override {}
    ZorderRTree(ZorderRTree const&) = default;
    ZorderRTree(ZorderRTree &&)     = default;
    ZorderRTree& operator=(ZorderRTree const&) = default;
    ZorderRTree& operator=(ZorderRTree &&)     = default;

    explicit ZorderRTree(const real_type margin)
        : cutoff_(0), margin_(margin), current_margin_(-1)
    {}

    bool valid() const noexcept override
    {
        return current_margin_ >= 0.0;
    }

    void initialize(neighbor_list_type& neighbors,
                    const system_type& sys, const potential_type& pot) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        if( ! is_simulator_traits<traits_type>::value)
        {
            MJOLNIR_LOG_WARN("spatial_partition.RTree is implemented for single"
                             " core simulation. Consider using CellList for"
                             " parallel simulation.");
        }

        const real_type max_cutoff = pot.max_cutoff_length();
        this->set_cutoff(max_cutoff);

        MJOLNIR_LOG_INFO(pot.name(), " cutoff = ", max_cutoff);

        this->make(neighbors, sys, pot);
        return;
    }

    void make(neighbor_list_type& neighbors,
              const system_type& sys, const potential_type& pot) override;

    bool reduce_margin(neighbor_list_type& neighbors, const real_type dmargin,
                       const system_type& sys, const potential_type& pot) override
    {
        this->current_margin_ -= dmargin;
        if(this->current_margin_ < 0)
        {
            this->make(neighbors, sys, pot);
            return true;
        }
        return false;
    }
    bool scale_margin(neighbor_list_type& neighbors, const real_type scale,
                const system_type& sys, const potential_type& pot) override
    {
        this->current_margin_ = (cutoff_ + current_margin_) * scale - cutoff_;
        if(this->current_margin_ < 0)
        {
            this->make(neighbors, sys, pot);
            return true;
        }
        return false;
    }

    real_type cutoff() const noexcept override {return this->cutoff_;}
    real_type margin() const noexcept override {return this->margin_;}

    base_type* clone() const override
    {
        return new ZorderRTree(margin_);
    }

  private:

    static std::uint32_t expand_bits(std::uint32_t v) noexcept
    {
        v = (v * 0x00010001u) & 0xFF0000FFu;
        v = (v * 0x00000101u) & 0x0F00F00Fu;
        v = (v * 0x00000011u) & 0xC30C30C3u;
        v = (v * 0x00000005u) & 0x49249249u;
        return v;
    }

    static std::uint32_t zindex(coordinate_type r) noexcept
    {
        math::X(r) = math::clamp<real_type>(math::X(r) * 1024, 0, 1023);
        math::Y(r) = math::clamp<real_type>(math::Y(r) * 1024, 0, 1023);
        math::Z(r) = math::clamp<real_type>(math::Z(r) * 1024, 0, 1023);
        const auto xx = expand_bits(static_cast<std::uint32_t>(math::X(r)));
        const auto yy = expand_bits(static_cast<std::uint32_t>(math::Y(r)));
        const auto zz = expand_bits(static_cast<std::uint32_t>(math::Z(r)));
        return xx * 4 + yy * 2 + zz;
    }

    static AABB merge_aabb(const AABB& box, const coordinate_type& p) noexcept
    {
        AABB merged(box);
        math::X(merged.lower) = std::min(math::X(merged.lower), math::X(p));
        math::Y(merged.lower) = std::min(math::Y(merged.lower), math::Y(p));
        math::Z(merged.lower) = std::min(math::Z(merged.lower), math::Z(p));

        math::X(merged.upper) = std::max(math::X(merged.upper), math::X(p));
        math::Y(merged.upper) = std::max(math::Y(merged.upper), math::Y(p));
        math::Z(merged.upper) = std::max(math::Z(merged.upper), math::Z(p));
        return merged;
    }
    static AABB merge_aabb(const AABB& lhs, const AABB& rhs) noexcept
    {
        AABB merged;
        math::X(merged.lower) = std::min(math::X(lhs.lower), math::X(rhs.lower));
        math::Y(merged.lower) = std::min(math::Y(lhs.lower), math::Y(rhs.lower));
        math::Z(merged.lower) = std::min(math::Z(lhs.lower), math::Z(rhs.lower));

        math::X(merged.upper) = std::max(math::X(lhs.upper), math::X(rhs.upper));
        math::Y(merged.upper) = std::max(math::Y(lhs.upper), math::Y(rhs.upper));
        math::Z(merged.upper) = std::max(math::Z(lhs.upper), math::Z(rhs.upper));
        return merged;
    }

    static std::pair<coordinate_type, coordinate_type> scaling(
        const system_type& sys, const std::vector<std::size_t>& participants,
        const UnlimitedBoundary<real_type, coordinate_type>&) noexcept
    {
        constexpr auto inf = std::numeric_limits<real_type>::infinity();
        AABB whole;
        whole.lower = math::make_coordinate<coordinate_type>( inf,  inf,  inf);
        whole.upper = math::make_coordinate<coordinate_type>(-inf, -inf, -inf);

        for(const auto& i : participants)
        {
            whole = merge_aabb(whole, sys.position(i));
        }
        const auto width = whole.upper - whole.lower;

        coordinate_type scale;
        math::X(scale) = real_type(1) / math::X(width);
        math::Y(scale) = real_type(1) / math::Y(width);
        math::Z(scale) = real_type(1) / math::Z(width);
        return std::make_pair(whole.lower, scale);
    }
    static std::pair<coordinate_type, coordinate_type> scaling(
        const system_type&, const std::vector<std::size_t>&,
        const CuboidalPeriodicBoundary<real_type, coordinate_type>& boundary) noexcept
    {
        coordinate_type scale;
        math::X(scale) = static_cast<real_type>(1) / math::X(boundary.width());
        math::Y(scale) = static_cast<real_type>(1) / math::Y(boundary.width());
        math::Z(scale) = static_cast<real_type>(1) / math::Z(boundary.width());
        return std::make_pair(boundary.lower_bound(), scale);
    }

    static bool overlaps(
        const coordinate_type& p, const real_type r, const AABB& box,
        const UnlimitedBoundary<real_type, coordinate_type>&) noexcept
    {
        const auto& l = box.lower;
        const auto& u = box.upper;
        if(math::X(p) + r < math::X(l) || math::X(u) < math::X(p) - r){return false;}
        if(math::Y(p) + r < math::Y(l) || math::Y(u) < math::Y(p) - r){return false;}
        if(math::Z(p) + r < math::Z(l) || math::Z(u) < math::Z(p) - r){return false;}
        return true;
    }
    static bool overlaps(
        const coordinate_type& p, const real_type r, const AABB& box,
        const CuboidalPeriodicBoundary<real_type, coordinate_type>& boundary) noexcept
    {
        const auto center = (box.upper + box.lower) * 0.5;
        const auto width  = (box.upper - box.lower) * 0.5;
        const auto dist = boundary.adjust_direction(center, p);

        if(math::X(width) + r < std::abs(math::X(dist))) {return false;}
        if(math::Y(width) + r < std::abs(math::Y(dist))) {return false;}
        if(math::Z(width) + r < std::abs(math::Z(dist))) {return false;}
        return true;
    }

    void set_cutoff(const real_type c) noexcept
    {
        this->cutoff_ = c;
        return;
    }
    void set_margin(const real_type m) noexcept
    {
        this->margin_ = m;
        return;
    }

    void diagnosis(const system_type& sys, const potential_type& pot) const;

  private:

    real_type cutoff_;
    real_type margin_;
    real_type current_margin_;
    std::vector<Node> tree_;
};

template<typename traitsT, typename PotentialT, std::size_t MaxElem>
void ZorderRTree<traitsT, PotentialT, MaxElem>::make(neighbor_list_type& neighbors,
        const system_type& sys, const potential_type& pot)
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();

    constexpr auto nil = std::numeric_limits<std::size_t>::max();

    // `participants` is a list that contains indices of particles that are
    // related to the potential.
    const auto& participants = pot.participants();

    neighbors.clear();
    tree_.clear();
    tree_.reserve(participants.size() * MaxElem / (MaxElem - 1));

    // ------------------------------------------------------------------------
    // calculate z-index

    const auto  scale  = scaling(sys, participants, sys.boundary());
    const auto& lower  = scale.first;
    const auto& rwidth = scale.second;

    // normally there are not so many particles
    assert(participants.size() < std::numeric_limits<std::uint32_t>::max());
    std::vector<std::pair<std::uint32_t, std::uint32_t>> zindices(participants.size());
    for(std::size_t i=0; i<participants.size(); ++i)
    {
        const auto idx = participants[i];
        auto pos = sys.position(idx) - lower;
        math::X(pos) *= math::X(rwidth);
        math::Y(pos) *= math::Y(rwidth);
        math::Z(pos) *= math::Z(rwidth);
        zindices[i] = std::make_pair(zindex(pos), static_cast<std::uint32_t>(idx));
    }
    // zindices contains {zindex, particle idx}

    // ------------------------------------------------------------------------
    // construct RTree

    std::sort(zindices.begin(), zindices.end());

    for(std::size_t i=0; i<zindices.size(); i += MaxElem)
    {
        Node node;
        node.parent  = nil;
        node.is_leaf = true;

        const auto idx = zindices[i].second;
        node.box.lower = sys.position(idx);
        node.box.upper = sys.position(idx);
        node.children.push_back(idx);

        for(std::size_t j=i+1; j<std::min(i + MaxElem, zindices.size()); ++j)
        {
            const auto jdx = zindices[j].second;
            node.box = merge_aabb(node.box, sys.position(jdx));
            node.children.push_back(jdx);
        }
        tree_.push_back(std::move(node));
    }

    std::size_t node_front = 0;
    while(tree_.size() - node_front > MaxElem)
    {
        const std::size_t N = tree_.size() - node_front;

        zindices.resize(N);
        for(std::size_t i=node_front; i<tree_.size(); ++i)
        {
            auto pos = (tree_[i].box.upper + tree_[i].box.lower) * 0.5 - lower;
            math::X(pos) *= math::X(rwidth);
            math::Y(pos) *= math::Y(rwidth);
            math::Z(pos) *= math::Z(rwidth);
            zindices[i - node_front] = std::make_pair(zindex(pos), i);
        }
        // zindices contains {zindex, node idx}
        std::sort(zindices.begin(), zindices.end());

        // update next node_front before adding a node of the current level
        node_front = tree_.size();

        for(std::size_t i=0; i<N; i += MaxElem)
        {
            const auto parent_idx = this->tree_.size();

            Node node;
            node.parent  = nil;
            node.is_leaf = false;

            // the first child (box does not need to be merged)
            {
                const auto child_node_idx = zindices[i].second;
                node.box = tree_[child_node_idx].box;
                node.children.push_back(child_node_idx);
                tree_[child_node_idx].parent = parent_idx;
            }
            for(std::size_t j=i+1; j<std::min(i + MaxElem, N); ++j)
            {
                const auto child_node_idx = zindices[j].second;
                node.box = merge_aabb(node.box, tree_[child_node_idx].box);
                node.children.push_back(child_node_idx);
                tree_[child_node_idx].parent = parent_idx;
            }
            tree_.push_back(std::move(node));
        }
    }
    const auto root = this->tree_.size();
    {
        Node node;
        node.parent  = nil;
        node.is_leaf = false;
        node.box = tree_[node_front].box;
        node.children.push_back(node_front);
        for(std::size_t i=node_front+1; i < tree_.size(); ++i)
        {
            node.box = merge_aabb(node.box, tree_[i].box);
            node.children.push_back(i);
            tree_[i].parent = root;
        }
        tree_.push_back(node);
    }

//     this->diagnosis(sys, pot);

    // ------------------------------------------------------------------------
    // construct NeighborList

    MJOLNIR_LOG_DEBUG("cell list is updated");

    const real_type r_c  = cutoff_ * (1. + margin_);
    const real_type r_c2 = r_c * r_c;

    const auto leading_participants = pot.leading_participants();

    std::vector<std::size_t> next_node;
    std::vector<neighbor_type> partner;
    for(std::size_t idx=0; idx<leading_participants.size(); ++idx)
    {
        next_node.clear();
        partner.clear();
        const auto   i = leading_participants[idx];
        const auto& ri = sys.position(i);

        next_node.push_back(root);
        while( ! next_node.empty())
        {
            const auto& node = tree_[next_node.back()];
            next_node.pop_back();

            if(node.is_leaf)
            {
                for(const auto& j : node.children)
                {
                    MJOLNIR_LOG_DEBUG("looking particle (", i, ", ", j, ")");
                    if( ! pot.has_interaction(i, j))
                    {
                        continue;
                    }
                    // here we don't need to search `participants` because
                    // cell list contains only participants. non-related
                    // particles are already filtered.

                    const auto& rj = sys.position(j);
                    if(math::length_sq(sys.adjust_direction(ri, rj)) < r_c2)
                    {
                        MJOLNIR_LOG_DEBUG("add index", j, "to verlet list", i);
                        partner.emplace_back(j, pot.prepare_params(i, j));
                    }
                }
            }
            else // internal node
            {
                for(const auto& child : node.children)
                {
                    if(overlaps(ri, r_c, tree_[child].box, sys.boundary()))
                    {
                        next_node.push_back(child);
                    }
                }
            }
        }

        // make the result consistent with NaivePairCalculation...
        std::sort(partner.begin(), partner.end());
        neighbors.add_list_for(i, partner.begin(), partner.end());

        // approximate the memory usage and avoid frequent memory allocation.
        // This block is just for efficiency, and does not affect on the result.
        if(idx == 16)
        {
            neighbors.reserve(neighbors.num_neighbors() / 16,
                              leading_participants.size(), sys.size());
        }
    }

    this->current_margin_ = cutoff_ * margin_;
    return ;
}

template<typename traitsT, typename PotentialT, std::size_t MaxElem>
void ZorderRTree<traitsT, PotentialT, MaxElem>::diagnosis(
        const system_type& sys, const potential_type& pot) const
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    std::vector<bool> found(sys.size(), false);
    for(const auto& node : tree_)
    {
        if(node.is_leaf)
        {
            for(const auto& child : node.children)
            {
                if(found.at(child))
                {
                    MJOLNIR_LOG_ERROR(child, "-th particle found twise!");
                }
                found.at(child) = true;

                const auto& pos = sys.position(child);

                bool included = true;
                if(math::X(pos) < math::X(node.box.lower) || math::X(node.box.upper) < math::X(pos)) {included = false;}
                if(math::Y(pos) < math::Y(node.box.lower) || math::Y(node.box.upper) < math::Y(pos)) {included = false;}
                if(math::Z(pos) < math::Z(node.box.lower) || math::Z(node.box.upper) < math::Z(pos)) {included = false;}

                if( ! included)
                {
                    MJOLNIR_LOG_ERROR(child, "-th particle sticks out of box!");
                }
            }
        }
        else
        {
            for(const auto& child : node.children)
            {
                const auto& child_box = tree_.at(child).box;

                bool included = true;
                if(math::X(child_box.lower) < math::X(node.box.lower) || math::X(node.box.upper) < math::X(child_box.upper)) {included = false;}
                if(math::Y(child_box.lower) < math::Y(node.box.lower) || math::Y(node.box.upper) < math::Y(child_box.upper)) {included = false;}
                if(math::Z(child_box.lower) < math::Z(node.box.lower) || math::Z(node.box.upper) < math::Z(child_box.upper)) {included = false;}

                if( ! included)
                {
                    MJOLNIR_LOG_ERROR(child, "-th node sticks out of box!");
                }
            }
        }
    }

    for(const auto& p : pot.participants())
    {
        if( ! found.at(p))
        {
            MJOLNIR_LOG_ERROR(p, "-th particle not found in the tree");
        }
    }
    return ;
}

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/forcefield/global/DebyeHuckelPotential.hpp>
#include <mjolnir/forcefield/global/ExcludedVolumePotential.hpp>
#include <mjolnir/forcefield/global/LennardJonesPotential.hpp>
#include <mjolnir/forcefield/global/UniformLennardJonesPotential.hpp>

namespace mjolnir
{
extern template class ZorderRTree<SimulatorTraits<double, CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ZorderRTree<SimulatorTraits<float,  CuboidalPeriodicBoundary>, DebyeHuckelPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ZorderRTree<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ZorderRTree<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ZorderRTree<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ZorderRTree<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ZorderRTree<SimulatorTraits<double, CuboidalPeriodicBoundary>, UniformLennardJonesPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class ZorderRTree<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UniformLennardJonesPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;

extern template class ZorderRTree<SimulatorTraits<double, UnlimitedBoundary>, DebyeHuckelPotential<SimulatorTraits<double, UnlimitedBoundary>>>;
extern template class ZorderRTree<SimulatorTraits<float,  UnlimitedBoundary>, DebyeHuckelPotential<SimulatorTraits<float,  UnlimitedBoundary>>>;

extern template class ZorderRTree<SimulatorTraits<double, UnlimitedBoundary>, ExcludedVolumePotential<SimulatorTraits<double, UnlimitedBoundary>>>;
extern template class ZorderRTree<SimulatorTraits<float,  UnlimitedBoundary>, ExcludedVolumePotential<SimulatorTraits<float,  UnlimitedBoundary>>>;

extern template class ZorderRTree<SimulatorTraits<double, UnlimitedBoundary>, LennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>>>;
extern template class ZorderRTree<SimulatorTraits<float,  UnlimitedBoundary>, LennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>>>;

extern template class ZorderRTree<SimulatorTraits<double, UnlimitedBoundary>, UniformLennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>>>;
extern template class ZorderRTree<SimulatorTraits<float,  UnlimitedBoundary>, UniformLennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>>>;
}
#endif // SEPARATE_BUILD

#endif /* MJOLNIR_CORE_ZORDER_R_TREE_HPP */
