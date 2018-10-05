#ifndef MJOLNIR_IGNORE_CHAIN_HPP
#define MJOLNIR_IGNORE_CHAIN_HPP
#include <mjolnir/util/throw_exception.hpp>
#include <string>
#include <memory>

// In some case, we need to ignore intra- or inter-chain interaction.
// For example, consider you have an elastic network model and charges on it
// that are tuned to reproduce the electrostatic potentials around the native
// structure. In that case, the electrostatic interaction should not be applied
// for the intra-chain pairs to avoid double-counting because the elastic
// network potential is modeled as a sum of all the intra-chain interactions.
//     IgnoreChain provides you a way to specify the rule to determine ignored
// pairs of particles. If you specify IgnoreSelf, all the intra-chain
// interactions would be ignored. IgnoreOthers ignores inter-chain interactions.
// Clearly, IgnoreNothing ignores nothing.

namespace mjolnir
{

template<typename ChainID>
struct IgnoreChainBase
{
    IgnoreChainBase()          = default;
    virtual ~IgnoreChainBase() = default;

    virtual bool is_ignored(const ChainID& i, const ChainID& j) const noexcept = 0;
    virtual const char* name() const noexcept = 0;
};

template<typename ChainID>
struct IgnoreNothing: public IgnoreChainBase<ChainID>
{
    IgnoreNothing()           = default;
    ~IgnoreNothing() override = default;

    bool is_ignored(const ChainID& i, const ChainID& j) const noexcept override
    {
        return false;
    }
    const char* name() const noexcept override {return "Nothing";}
};

template<typename ChainID>
struct IgnoreSelf: public IgnoreChainBase<ChainID>
{
    IgnoreSelf()           = default;
    ~IgnoreSelf() override = default;

    bool is_ignored(const ChainID& i, const ChainID& j) const noexcept override
    {
        return i == j;
    }
    const char* name() const noexcept override {return "Self";}
};

template<typename ChainID>
struct IgnoreOthers: public IgnoreChainBase<ChainID>
{
    IgnoreOthers()           = default;
    ~IgnoreOthers() override = default;

    bool is_ignored(const ChainID& i, const ChainID& j) const noexcept override
    {
        return i != j;
    }
    const char* name() const noexcept override {return "Others";}
};

template<typename ChainID>
struct IgnoreChain
{
    IgnoreChain(const std::string& name): ignore_chain_(nullptr)
    {
        this->reset(name);
    }

    template<typename IgnoreSomething>
    IgnoreChain(std::unique_ptr<IgnoreSomething>&& ptr)
        : ignore_chain_(std::move(ptr))
    {}
    ~IgnoreChain() = default;

    IgnoreChain(IgnoreChain const& rhs) {this->reset(rhs.name());}
    IgnoreChain(IgnoreChain&&      rhs) = default;
    IgnoreChain& operator=(IgnoreChain const& rhs) {this->reset(rhs.name()); return *this;}
    IgnoreChain& operator=(IgnoreChain&&      rhs) = default;

    bool is_ignored(const ChainID& i, const ChainID& j) const noexcept
    {
        return ignore_chain_->is_ignored(i, j);
    }
    const char* name() const noexcept {return ignore_chain_->name();}

    void reset(const std::string& name)
    {
        if(name == "Nothing")
        {
            ignore_chain_.reset(new IgnoreNothing<ChainID>);
        }
        else if(name == "Self")
        {
            ignore_chain_.reset(new IgnoreSelf<ChainID>);
        }
        else if(name == "Others")
        {
            ignore_chain_.reset(new IgnoreOthers<ChainID>);
        }
        else
        {
            throw_exception<std::invalid_argument>("IgnoreChain::IgnoreChain: ",
                "unknown signeture appreaed: `" + name +
                "`. Only `Nothing`, `Self`, or `Others` are allowed.");
        }
    }

  private:

    std::unique_ptr<IgnoreChainBase<ChainID>> ignore_chain_;
};

} // mjolnir
#endif// MJOLNIR_IGNORE_CHAIN_HPP
