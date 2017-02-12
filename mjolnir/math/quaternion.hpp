#ifndef MJOLNIR_QUATERNION
#define MJOLNIR_QUATERNION
#include <type_traits>
#include <array>
#include <cmath>

namespace mjolnir
{

template<typename realT>
class quaternion
{
  public:
    typedef realT real_type;
    typedef real_type value_type;
    typedef std::array<real_type, 4> container_type;

  public:
    quaternion()  = default;
    ~quaternion() = default;
    quaternion(const quaternion&) = default;
    quaternion(quaternion&&)      = default;
    quaternion& operator=(const quaternion&) = default;
    quaternion& operator=(quaternion&&)      = default;

    template<typename T1, typename T2, typename T3, typename T4, class =
        typename std::enable_if<std::is_convertible<T1, real_type>::value &&
            std::is_convertible<T2, real_type>::value &&
            std::is_convertible<T3, real_type>::value &&
            std::is_convertible<T4, real_type>::value>::type>
    quaternion(T1&& v1, T2&& v2, T3&& v3, T4&& v4)
    : values_{{std::forward<real_type>(v1), std::forward<real_type>(v2),
               std::forward<real_type>(v3), std::forward<real_type>(v4)}}
    {}

    template<typename T1, typename T2, typename T3, typename T4, class =
        typename std::enable_if<std::is_convertible<T1, real_type>::value &&
            std::is_convertible<T2, real_type>::value &&
            std::is_convertible<T3, real_type>::value &&
            std::is_convertible<T4, real_type>::value>::type>
    quaternion(const T1& v1, const T2& v2, const T3& v3, const T4& v4)
    : values_{{static_cast<real_type>(v1), static_cast<real_type>(v2),
               static_cast<real_type>(v3), static_cast<real_type>(v4)}}
    {}

    template<typename T, typename std::enable_if<
        std::is_convertible<T, real_type>::value, std::nullptr_t>::type = nullptr>
    quaternion& operator+=(const quaternion<T>&) noexcept;
    template<typename T, typename std::enable_if<
        std::is_convertible<T, real_type>::value, std::nullptr_t>::type = nullptr>
    quaternion& operator-=(const quaternion<T>&) noexcept;
    template<typename T, typename std::enable_if<
        std::is_convertible<T, real_type>::value, std::nullptr_t>::type = nullptr>
    quaternion& operator*=(const quaternion<T>&) noexcept;
    template<typename T, typename std::enable_if<
        std::is_convertible<T, real_type>::value, std::nullptr_t>::type = nullptr>
    quaternion& operator/=(const quaternion<T>&) noexcept;
    template<typename T, typename std::enable_if<
        std::is_convertible<T, real_type>::value, std::nullptr_t>::type = nullptr>
    quaternion& operator*=(const T&) noexcept;
    template<typename T, typename std::enable_if<
        std::is_convertible<T, real_type>::value, std::nullptr_t>::type = nullptr>
    quaternion& operator/=(const T&) noexcept;

    real_type& operator[](const std::size_t i)       {return values_[i];}
    real_type  operator[](const std::size_t i) const {return values_[i];}
    real_type& at(const std::size_t i)       {return values_.at(i);}
    real_type  at(const std::size_t i) const {return values_.at(i);}

  private:

    container_type values_;
};

template<typename realT>
inline realT norm_sq(const quaternion<realT>& q)
{
    return q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
}

template<typename realT>
inline realT norm(const quaternion<realT>& q)
{
    return std::sqrt(norm_sq(q));
}

template<typename realT>
inline quaternion<realT> conj(const quaternion<realT>& q)
{
    return quaternion<realT>(q[0], -q[1], -q[2], -q[3]);
}

template<typename realT>
template<typename T, typename std::enable_if<
    std::is_convertible<T, realT>::value, std::nullptr_t>::type>
inline quaternion<realT>&
quaternion<realT>::operator+=(const quaternion<T>& rhs) noexcept
{
    values_[0] += rhs[0];
    values_[1] += rhs[1];
    values_[2] += rhs[2];
    values_[3] += rhs[3];
    return *this;
}

template<typename realT>
template<typename T, typename std::enable_if<
    std::is_convertible<T, realT>::value, std::nullptr_t>::type>
inline quaternion<realT>&
quaternion<realT>::operator-=(const quaternion<T>& rhs) noexcept
{
    values_[0] -= rhs[0];
    values_[1] -= rhs[1];
    values_[2] -= rhs[2];
    values_[3] -= rhs[3];
    return *this;
}

template<typename realT>
template<typename T, typename std::enable_if<
    std::is_convertible<T, realT>::value, std::nullptr_t>::type>
inline quaternion<realT>&
quaternion<realT>::operator*=(const quaternion<T>& rhs) noexcept
{
    *this = (*this) * rhs;
    return *this;
}

template<typename realT>
template<typename T, typename std::enable_if<
    std::is_convertible<T, realT>::value, std::nullptr_t>::type>
inline quaternion<realT>&
quaternion<realT>::operator/=(const quaternion<T>& rhs) noexcept
{
    *this = (*this) / rhs;
    return *this;
}

template<typename realT>
template<typename T, typename std::enable_if<
    std::is_convertible<T, realT>::value, std::nullptr_t>::type>
inline quaternion<realT>&
quaternion<realT>::operator*=(const T& rhs) noexcept
{
    values_[0] *= rhs;
    values_[1] *= rhs;
    values_[2] *= rhs;
    values_[3] *= rhs;
    return *this;
}

template<typename realT>
template<typename T, typename std::enable_if<
    std::is_convertible<T, realT>::value, std::nullptr_t>::type>
inline quaternion<realT>&
quaternion<realT>::operator/=(const T& rhs) noexcept
{
    values_[0] /= rhs;
    values_[1] /= rhs;
    values_[2] /= rhs;
    values_[3] /= rhs;
    return *this;
}

template<typename T1, typename T2, typename std::enable_if<
    std::is_convertible<T1, T2>::value, std::nullptr_t>::type = nullptr>
inline quaternion<decltype(std::declval<T1>() + std::declval<T2>())>
operator+(const quaternion<T1>& lhs, const quaternion<T2>& rhs)
{
    return quaternion<decltype(std::declval<T1>() + std::declval<T2>())>(
            lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2], lhs[3] + rhs[3]);
}

template<typename T1, typename T2, typename std::enable_if<
    std::is_convertible<T1, T2>::value, std::nullptr_t>::type = nullptr>
inline quaternion<decltype(std::declval<T1>() - std::declval<T2>())>
operator-(const quaternion<T1>& lhs, const quaternion<T2>& rhs)
{
    return quaternion<decltype(std::declval<T1>() - std::declval<T2>())>(
            lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2], lhs[3] - rhs[3]);
}

template<typename T1, typename T2, typename std::enable_if<
    std::is_convertible<T1, T2>::value, std::nullptr_t>::type = nullptr>
inline quaternion<decltype(std::declval<T1>() * std::declval<T2>())>
operator*(const quaternion<T1>& lhs, const quaternion<T2>& rhs)
{
    typedef decltype(std::declval<T1>() * std::declval<T2>()) realT;
    const realT a = lhs[0]*rhs[0] - lhs[1]*rhs[1] - lhs[2]*rhs[2] - lhs[3]*rhs[3];
    const realT b = lhs[0]*rhs[1] + lhs[1]*rhs[0] + lhs[2]*rhs[3] - lhs[3]*rhs[2];
    const realT c = lhs[0]*rhs[2] - lhs[1]*rhs[3] + lhs[2]*rhs[0] + lhs[3]*rhs[1];
    const realT d = lhs[0]*rhs[3] + lhs[1]*rhs[2] - lhs[2]*rhs[1] + lhs[3]*rhs[0];
    return quaternion<realT>(a, b, c, d);
}

template<typename T1, typename T2, typename std::enable_if<
    std::is_convertible<T1, T2>::value, std::nullptr_t>::type = nullptr>
inline quaternion<decltype(std::declval<T1>() / std::declval<T2>())>
operator/(const quaternion<T1>& lhs, const quaternion<T2>& rhs)
{
    typedef decltype(std::declval<T1>() / std::declval<T2>()) realT;
    return lhs * conj(rhs) / norm_sq(rhs);
}

template<typename T1, typename T2, typename std::enable_if<
    std::is_convertible<T1, T2>::value, std::nullptr_t>::type = nullptr>
quaternion<decltype(std::declval<T1>() * std::declval<T2>())>
operator*(const quaternion<T1>& lhs, const T2& rhs)
{
    return quaternion<decltype(std::declval<T1>() * std::declval<T2>())>(
            lhs[0] * rhs, lhs[1] * rhs, lhs[2] * rhs, lhs[3] * rhs);
}

template<typename T1, typename T2, typename std::enable_if<
    std::is_convertible<T1, T2>::value, std::nullptr_t>::type = nullptr>
quaternion<decltype(std::declval<T1>() * std::declval<T2>())>
operator*(const T1& lhs, const quaternion<T2>& rhs)
{
    return quaternion<decltype(std::declval<T1>() * std::declval<T2>())>(
            lhs * rhs[0], lhs * rhs[1], lhs * rhs[2], lhs * rhs[3]);
}

template<typename T1, typename T2, typename std::enable_if<
    std::is_convertible<T1, T2>::value, std::nullptr_t>::type = nullptr>
quaternion<decltype(std::declval<T1>() / std::declval<T2>())>
operator/(const quaternion<T1>& lhs, const T2& rhs)
{
    return quaternion<decltype(std::declval<T1>() / std::declval<T2>())>(
            lhs[0] / rhs, lhs[1] / rhs, lhs[2] / rhs, lhs[3] / rhs);
}

} // mjolnir
#endif /*MJOLNIR_QUATERNION*/
