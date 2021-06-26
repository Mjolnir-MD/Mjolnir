#ifndef MJOLNIR_CORE_DYNAMIC_VARIABLE_HPP
#define MJOLNIR_CORE_DYNAMIC_VARIABLE_HPP
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/is_finite.hpp>
#include <limits>
#include <memory>

//
// Dynamic variable is a (non-) physical parameter in a system that has
// (virtual) mass, velocity and force and updated by integrator.
//
namespace mjolnir
{

template<typename realT>
struct DynamicVariableBase
{
    using real_type = realT;

    real_type x, v, f;
    real_type m, gamma;
    real_type lower, upper;

    DynamicVariableBase(real_type x_, real_type v_, real_type f_,
                        real_type m_, real_type gamma_,
                        real_type lower_, real_type upper_)
        : x(x_), v(v_), f(f_), m(m_), gamma(gamma_), lower(lower_), upper(upper_)
    {}
    virtual ~DynamicVariableBase() = default;
    virtual void update(real_type x, real_type v, real_type f) noexcept = 0;
    virtual DynamicVariableBase<real_type>* clone() const = 0;
};

// It makes `DynamicVariable`s copyable and easier to handle.
template<typename realT>
struct DynamicVariable
{
    using real_type = typename DynamicVariableBase<realT>::real_type;

    template<typename DynVar>
    explicit DynamicVariable(DynVar&& dynvar)
        : resource_(make_unique<DynamicVariableBase<real_type>>(std::forward<DynVar>(dynvar)))
    {
        static_assert(std::is_base_of<DynamicVariableBase<realT>, DynVar>::value, "");
    }
    DynamicVariable(): resource_(nullptr) {}

    DynamicVariable(const DynamicVariable& other)
        : resource_(other.resource_->clone())
    {}
    DynamicVariable& operator=(const DynamicVariable& other)
    {
        this->resource_.reset(other.resource_->clone());
        return *this;
    }
    DynamicVariable(DynamicVariable&&) = default;
    DynamicVariable& operator=(DynamicVariable&&) = default;
    ~DynamicVariable() = default;

    real_type x() const noexcept {return resource_->x;}
    real_type v() const noexcept {return resource_->v;}
    real_type f() const noexcept {return resource_->f;}
    real_type m() const noexcept {return resource_->m;}
    real_type gamma() const noexcept {return resource_->gamma;}

    void update(real_type x, real_type v, real_type f) const noexcept
    {
        resource_->update(x, v, f);
    }

  private:
    std::unique_ptr<DynamicVariableBase<realT>> resource_;
};

template<typename realT>
struct DefaultDynamicVariable : public DynamicVariableBase<realT>
{
    using base_type = DynamicVariableBase<realT>;
    using real_type = typename base_type::real_type;

    DefaultDynamicVariable(const real_type x, const real_type v, const real_type f,
                           const real_type m, const real_type gamma)
        : base_type{x, v, f, m, gamma,
                    -std::numeric_limits<real_type>::infinity(),
                     std::numeric_limits<real_type>::infinity()}
    {}
    ~DefaultDynamicVariable() override = default;

    void update(real_type x, real_type v, real_type f) noexcept override
    {
        this->x = x;
        this->v = v;
        this->f = f;
        return;
    }
    DynamicVariableBase<real_type>* clone() const override
    {
        return new DefaultDynamicVariable(*this);
    }
};

template<typename realT>
struct PeriodicDynamicVariable : public DynamicVariableBase<realT>
{
    using base_type = DynamicVariableBase<realT>;
    using real_type = typename base_type::real_type;

    PeriodicDynamicVariable(const real_type x, const real_type v, const real_type f,
                            const real_type m, const real_type gamma,
                            const real_type lower, const real_type upper)
        : base_type(x, v, f, m, gamma, lower, upper)
    {}
    ~PeriodicDynamicVariable() override = default;

    void update(real_type x, real_type v, real_type f) noexcept override
    {
        this->x = x;
        if      (this->x <  this->lower) {this->x += (this->upper - this->lower);}
        else if (this->upper <= this->x) {this->x -= (this->upper - this->lower);}
        this->v = v;
        this->f = f;
        return;
    }
    DynamicVariableBase<real_type>* clone() const override
    {
        return new PeriodicDynamicVariable(*this);
    }
};

template<typename realT>
struct RepulsiveDynamicVariable : public DynamicVariableBase<realT>
{
    using base_type = DynamicVariableBase<realT>;
    using real_type = typename base_type::real_type;

    RepulsiveDynamicVariable(const real_type x, const real_type v, const real_type f,
                             const real_type m, const real_type gamma,
                             const real_type lower, const real_type upper)
        : base_type(x, v, f, m, gamma, lower, upper)
    {}
    ~RepulsiveDynamicVariable() override = default;

    void update(real_type x, real_type v, real_type f) noexcept override
    {
        if(x < this->lower)
        {
            this->x = 2 * this->lower - x; // lower + (lower - x)
            this->v = -v;
        }
        else if (this->upper < x)
        {
            this->x = 2 * this->upper - x; // upper - (x - upper)
            this->v = -v;
        }
        else
        {
            this->x = x;
            this->v = v;
        }
        this->f = f;
        return;
    }
    DynamicVariableBase<real_type>* clone() const override
    {
        return new RepulsiveDynamicVariable(*this);
    }
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class DynamicVariableBase<double>;
extern template class DynamicVariableBase<float >;
extern template class DynamicVariable<double>;
extern template class DynamicVariable<float >;
extern template class DefaultDynamicVariable<double>;
extern template class DefaultDynamicVariable<float >;
extern template class PeriodicDynamicVariable<double>;
extern template class PeriodicDynamicVariable<float >;
extern template class RepulsiveDynamicVariable<double>;
extern template class RepulsiveDynamicVariable<float >;
#endif

} // mjolnir
#endif// MJOLNIR_CORE_DYNAMIC_VARIABLE_HPP
