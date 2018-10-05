#ifndef MJOLNIR_CORE_UNIT_HPP
#define MJOLNIR_CORE_UNIT_HPP
#include <ratio>

namespace mjolnir
{

namespace unit
{

// -----------------------------------------------------------------------------
// generalized ratio class to use 1/mol or something that cannot be represented
// by a rational in int64_t (like `per mol`). The result of this might include
// numerical errors because it uses floating-point arithmetics!
// -----------------------------------------------------------------------------

template<typename Numer, typename Denom>
struct Ratio
{
    template<typename T>
    struct value_type_of {using type = typename T::value_type;};
    template<std::intmax_t N1, std::intmax_t D1>
    struct value_type_of<std::ratio<N1, D1>> {using type = double;};
    template<typename T>
    using value_type_of_t = typename value_type_of<T>::type;

    template<typename T>
    struct value_of {static constexpr value_type_of_t<T> value = T::value;};
    template<std::intmax_t N1, std::intmax_t D1>
    struct value_of<std::ratio<N1, D1>>
    {
        static constexpr double value = static_cast<double>(N1) / D1;
    };

    using numerator_value_type   = value_type_of_t<Numer>;
    using denominator_value_type = value_type_of_t<Denom>;

    static constexpr numerator_value_type   num = value_of<Numer>::value;
    static constexpr denominator_value_type den = value_of<Denom>::value;

    using value_type = decltype(num / den);
    static constexpr value_type value = num / den;
};
template<typename Numer, typename Denom>
constexpr typename Ratio<Numer, Denom>::numerator_value_type Ratio<Numer, Denom>::num;
template<typename Numer, typename Denom>
constexpr typename Ratio<Numer, Denom>::denominator_value_type Ratio<Numer, Denom>::den;
template<typename Numer, typename Denom>
constexpr typename Ratio<Numer, Denom>::value_type Ratio<Numer, Denom>::value;

// for ease
using One = std::integral_constant<std::intmax_t, 1>;

template<typename Real, typename FromRatio, typename ToRatio>
struct ratio_between
{
    using real_type = Real;

    static constexpr real_type value =
        (static_cast<real_type>(FromRatio::num) * static_cast<real_type>(ToRatio::den)) /
        (static_cast<real_type>(FromRatio::den) * static_cast<real_type>(ToRatio::num));
};
template<typename Real, typename FromRatio, typename ToRatio>
constexpr Real ratio_between<Real, FromRatio, ToRatio>::value;

// specialization for std::ratio<N, D>

template<typename Real, std::intmax_t N1, std::intmax_t D1,
                        std::intmax_t N2, std::intmax_t D2>
struct ratio_between<Real, std::ratio<N1, D1>, std::ratio<N2, D2>>
{
    using real_type  = Real;
    using ratio_type = std::ratio_divide<std::ratio<N1, D1>, std::ratio<N2, D2>>;

    static constexpr real_type value =
        static_cast<real_type>(ratio_type::num) /
        static_cast<real_type>(ratio_type::den);
};
template<typename Real, std::intmax_t N1, std::intmax_t D1,
                        std::intmax_t N2, std::intmax_t D2>
constexpr Real ratio_between<Real, std::ratio<N1,D1>, std::ratio<N2,D2>>::value;

// -------------------------------------------------------------------------
// conversion between Units (From -> To)
// -------------------------------------------------------------------------

template<typename Real, typename FromUnit, typename ToUnit>
struct convert
{
    using real_type = Real;

    static constexpr real_type coef = ratio_between<Real,
        typename FromUnit::ratio_type, typename ToUnit::ratio_type>::value;

    static constexpr real_type invoke(real_type from) noexcept
    {
        return from * coef;
    }
};
template<typename Real, typename FromUnit, typename ToUnit>
constexpr Real convert<Real, FromUnit, ToUnit>::coef;

// -------------------------------------------------------------------------
// tag (+ratio) types to represent physical units
// -------------------------------------------------------------------------

template<typename T, typename R = std::ratio<1, 1>>
struct meter
{
    using value_type = T;
    using ratio_type = R;
};

template<typename T> using micro_meter = meter<T, std::micro>;
template<typename T> using nano_meter  = meter<T, std::nano>;
template<typename T> using angstrom    = meter<T, std::ratio<1, 10000000000>>;
template<typename T> using pico_meter  = meter<T, std::pico>;

template<typename T, typename R = std::ratio<1, 1>>
struct mol
{
    using value_type = T;
    using ratio_type = R;

    // XXX mol can be used as a ratio later...
    static constexpr value_type value = static_cast<value_type>(6.022140857e23);
};
template<typename T, typename R>
constexpr T mol<T, R>::value;

template<typename T> using per_mol = Ratio<One, mol<T>>;

template<typename T, typename R = std::ratio<1, 1>>
struct joule
{
    using value_type = T;
    using ratio_type = R;
};


template<typename T> using kilo_joule         = joule<T, std::kilo>;
template<typename T> using joule_per_mol      = joule<T, per_mol<T>>;
template<typename T> using kilo_joule_per_mol = joule<T, Ratio<std::kilo, mol<T>>>;

template<typename T, typename R> using J          = joule<T, R>;
template<typename T>             using kJ         = kilo_joule<T>;
template<typename T>             using J_per_mol  = joule_per_mol<T>;
template<typename T>             using kJ_per_mol = kilo_joule_per_mol<T>;

template<typename T, typename R = std::ratio<1, 1>>
struct calorie
{
    using value_type = T;
    using ratio_type = R;
};

template<typename T> using kilo_calorie         = calorie<T, std::kilo>;
template<typename T> using calorie_per_mol      = calorie<T, per_mol<T>>;
template<typename T> using kilo_calorie_per_mol = calorie<T, Ratio<std::kilo, mol<T>>>;

template<typename T, typename R> using cal          = calorie<T, R>;
template<typename T>             using kcal         = kilo_calorie<T>;
template<typename T>             using cal_per_mol  = calorie_per_mol<T>;
template<typename T>             using kcal_per_mol = kilo_calorie_per_mol<T>;





template<typename T, template<typename> class Unit>
struct quantity
{
    using value_type = T;
    using unit_type  = Unit<value_type>;

    quantity(value_type v): value(v) {}

    quantity(): value(0){};
    ~quantity() = default;
    quantity(quantity const&) = default;
    quantity(quantity &&)     = default;
    quantity& operator=(quantity const&) = default;
    quantity& operator=(quantity &&)     = default;

    template<typename T2, template<typename> class Unit2>
    quantity(const quantity<T2, Unit2>& rhs)
        : value(convert<T, Unit2<value_type>, unit_type>::invoke(rhs.value))
    {}
    template<typename T2, template<typename> class Unit2>
    quantity& operator=(const quantity<T2, Unit2>& rhs) noexcept
    {
        this->value = convert<T, Unit2<value_type>, unit_type>::invoke(rhs.value);
        return *this;
    }

    template<typename T2, template<typename> class Unit2>
    quantity& operator+=(const quantity<T2, Unit2>& rhs) noexcept
    {
        this->value += convert<T, Unit2<value_type>, unit_type>::invoke(rhs.value);
        return *this;
    }

    value_type value;
};

} // unit
} // mjolnir
#endif// MJOLNIR_CORE_UNIT_HPP
