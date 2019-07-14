#ifndef MJOLNIR_POTENTIAL_GLOBAL_IGNORE_MOLECULE_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_IGNORE_MOLECULE_HPP
#include <mjolnir/util/throw_exception.hpp>
#include <string>
#include <memory>

// In some case, we need to ignore intra- or inter-molecule interaction.
// For example, consider you have an elastic network model and charges on it
// that are tuned to reproduce the electrostatic potentials around the native
// structure. In that case, the electrostatic interaction should not be applied
// for the intra-molecule pairs to avoid double-counting because the elastic
// network potential is modeled as a sum of all the intra-molecule interactions.
//     IgnoreMolecule provides you a way to specify the rule to determine ignored
// pairs of particles. If you specify IgnoreSelf, all the intra-molecule
// interactions would be ignored. IgnoreOthers ignores inter-molecule interactions.
// Clearly, IgnoreNothing ignores nothing.

namespace mjolnir
{

template<typename MoleculeID>
struct IgnoreMoleculeBase
{
    IgnoreMoleculeBase()          = default;
    virtual ~IgnoreMoleculeBase() = default;

    virtual bool is_ignored(const MoleculeID& i, const MoleculeID& j) const noexcept = 0;
    virtual const char* name() const noexcept = 0;
};

template<typename MoleculeID>
struct IgnoreNothing: public IgnoreMoleculeBase<MoleculeID>
{
    IgnoreNothing()           = default;
    ~IgnoreNothing() override = default;

    bool is_ignored(const MoleculeID&, const MoleculeID&) const noexcept override
    {
        return false;
    }
    const char* name() const noexcept override {return "Nothing";}
};

template<typename MoleculeID>
struct IgnoreSelf: public IgnoreMoleculeBase<MoleculeID>
{
    IgnoreSelf()           = default;
    ~IgnoreSelf() override = default;

    bool is_ignored(const MoleculeID& i, const MoleculeID& j) const noexcept override
    {
        return i == j;
    }
    const char* name() const noexcept override {return "Self";}
};

template<typename MoleculeID>
struct IgnoreOthers: public IgnoreMoleculeBase<MoleculeID>
{
    IgnoreOthers()           = default;
    ~IgnoreOthers() override = default;

    bool is_ignored(const MoleculeID& i, const MoleculeID& j) const noexcept override
    {
        return i != j;
    }
    const char* name() const noexcept override {return "Others";}
};

template<typename MoleculeID>
struct IgnoreMolecule
{
    explicit IgnoreMolecule(const std::string& name): ignore_mol_(nullptr)
    {
        this->reset(name);
    }

    template<typename IgnoreSomething>
    explicit IgnoreMolecule(std::unique_ptr<IgnoreSomething>&& ptr)
        : ignore_mol_(std::move(ptr))
    {}
    ~IgnoreMolecule() = default;

    IgnoreMolecule(IgnoreMolecule const& rhs) {this->reset(rhs.name());}
    IgnoreMolecule(IgnoreMolecule&&      rhs) = default;
    IgnoreMolecule& operator=(IgnoreMolecule const& rhs) {this->reset(rhs.name()); return *this;}
    IgnoreMolecule& operator=(IgnoreMolecule&&      rhs) = default;

    bool is_ignored(const MoleculeID& i, const MoleculeID& j) const noexcept
    {
        return ignore_mol_->is_ignored(i, j);
    }
    const char* name() const noexcept {return ignore_mol_->name();}

    void reset(const std::string& name)
    {
        if(name == "Nothing")
        {
            ignore_mol_.reset(new IgnoreNothing<MoleculeID>);
        }
        else if(name == "Self" || name == "Intra")
        {
            ignore_mol_.reset(new IgnoreSelf<MoleculeID>);
        }
        else if(name == "Others" || name == "Inter")
        {
            ignore_mol_.reset(new IgnoreOthers<MoleculeID>);
        }
        else
        {
            throw_exception<std::invalid_argument>("IgnoreMolecule::IgnoreMolecule: ",
                "unknown signeture appreaed: `" + name +
                "`. Only `Nothing`, `Self`, or `Others` are allowed.");
        }
    }

  private:

    std::unique_ptr<IgnoreMoleculeBase<MoleculeID>> ignore_mol_;
};

} // mjolnir
#endif// MJOLNIR_IGNORE_MOLECULE_HPP
