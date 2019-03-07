#ifndef MJOLNIR_MATH_MATRIX
#define MJOLNIR_MATH_MATRIX
#include <mjolnir/util/type_traits.hpp>
#include <algorithm>
#include <array>
#include <ostream>
#include <cmath>

namespace mjolnir
{
namespace math
{

// matrix type.
//    1   ...   C
//  1 x00 ... x0M    N = R-1
//  . x10 ... x1M    M = C-1
//  .   . ...   .
//  R xN0 ... xNM    R * C matrix
//
// Access by using `operator()` in the following way.
// Matrix<double, 3, 3> m;
// assert(m(1,0) == x10);
// assert(m(2,1) == x21);
//
// Note that the order of `operator[]` is not guaranteed.
// *DO NOT* use `mat[i]` to modify the element.
// It is only for the internal use.
template<typename realT, std::size_t R, std::size_t C>
class Matrix
{
  public:

    static constexpr std::size_t    row_size = R;
    static constexpr std::size_t column_size = C;
    static constexpr std::size_t  total_size = R * C;

    using value_type      = realT;
    using storage_type    = std::array<value_type, total_size>;
    using pointer         = value_type*;
    using const_pointer   = value_type const*;
    using reference       = value_type&;
    using const_reference = value_type const&;
    using size_type       = std::size_t;

  public:
    Matrix() : values_{{}}{}
    ~Matrix() = default;
    Matrix(Matrix const&) = default;
    Matrix(Matrix &&)     = default;
    Matrix& operator=(Matrix const&) = default;
    Matrix& operator=(Matrix &&)     = default;

    template<typename ... Ts, typename std::enable_if<
        sizeof...(Ts) == total_size, std::nullptr_t>::type = nullptr>
    Matrix(Ts&& ... vs) noexcept: values_{{static_cast<value_type>(vs)...}}
    {
        static_assert(sizeof...(Ts) == total_size, "");
        static_assert(conjunction<
                std::is_convertible<Ts, value_type> ...>::value, "");
    }

    template<typename T>
    Matrix(const std::array<T, total_size>& rhs) noexcept
    {
        static_assert(std::is_convertible<T, value_type>::value, "");
        for(std::size_t i=0; i<total_size; ++i){this->values_[i] = rhs[i];}
    }
    template<typename T>
    Matrix(const Matrix<T, R, C>& rhs) noexcept
    {
        static_assert(std::is_convertible<T, value_type>::value, "");
        for(std::size_t i=0; i<total_size; ++i){this->values_[i] = rhs[i];}
    }

    template<typename T>
    Matrix& operator=(const Matrix<T, R, C>& rhs) noexcept
    {
        static_assert(std::is_convertible<T, value_type>::value, "");
        for(std::size_t i=0; i<total_size; ++i){this->values_[i] = rhs[i];}
    }

    template<typename T>
    Matrix& operator+=(const Matrix<T, R, C>& rhs) noexcept
    {
        static_assert(std::is_convertible<T, value_type>::value, "");
        for(std::size_t i=0; i<total_size; ++i){this->values_[i] += rhs[i];}
        return *this;
    }
    template<typename T>
    Matrix& operator-=(const Matrix<T, R, C>& rhs) noexcept
    {
        static_assert(std::is_convertible<T, value_type>::value, "");
        for(std::size_t i=0; i<total_size; ++i){this->values_[i] -= rhs[i];}
        return *this;
    }

    template<typename T>
    Matrix& operator*=(const T& rhs) noexcept
    {
        static_assert(std::is_convertible<T, value_type>::value, "");
        for(std::size_t i=0; i<total_size; ++i){this->values_[i] *= rhs;}
        return *this;
    }
    template<typename T>
    Matrix& operator/=(const T& rhs) noexcept
    {
        static_assert(std::is_convertible<T, value_type>::value, "");
        for(std::size_t i=0; i<total_size; ++i){this->values_[i] /= rhs;}
        return *this;
    }

    value_type  at(size_type i, size_type j) const {return values_.at(i*C+j);}
    value_type& at(size_type i, size_type j)       {return values_.at(i*C+j);}
    value_type  operator()(size_type i, size_type j) const noexcept {return values_[i*C+j];}
    value_type& operator()(size_type i, size_type j)       noexcept {return values_[i*C+j];}

    value_type  at(size_type i) const {return values_.at(i);}
    value_type& at(size_type i)       {return values_.at(i);}
    value_type  operator[](size_type i) const noexcept {return values_[i];}
    value_type& operator[](size_type i)       noexcept {return values_[i];}

    pointer       data()       noexcept {return values_.data();}
    const_pointer data() const noexcept {return values_.data();}

    bool diagnosis() const noexcept {return true;}

    void zero() noexcept {values_.fill(value_type(0));}

  private:
    storage_type values_;
};
template<typename realT, std::size_t R, std::size_t C>
constexpr std::size_t Matrix<realT, R, C>::row_size;
template<typename realT, std::size_t R, std::size_t C>
constexpr std::size_t Matrix<realT, R, C>::column_size;
template<typename realT, std::size_t R, std::size_t C>
constexpr std::size_t Matrix<realT, R, C>::total_size;

template<typename charT, typename traitsT,
         typename T, std::size_t R, std::size_t C>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const Matrix<T, R, C>& mat)
{
    for(std::size_t i=0; i<R; ++i)
    {
        os << '(';
        for(std::size_t j=0; j<R; ++j)
        {
            if(j!=0) {os << ',';}
            os << mat(i, j);
        }
        os << ")\n";
    }
    return os;
}

// negation operator
template<typename T, std::size_t R, std::size_t C>
inline Matrix<T, R, C>
operator-(const Matrix<T, R, C>& rhs) noexcept
{
    Matrix<T, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i) {retval[i] = -rhs[i];}
    return retval;
}

template<typename T, std::size_t R, std::size_t C>
inline Matrix<T, R, C>
operator+(const Matrix<T, R, C>& lhs, const Matrix<T, R, C>& rhs) noexcept
{
    Matrix<T, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i) {retval[i] = lhs[i] + rhs[i];}
    return retval;
}
template<typename T, std::size_t R, std::size_t C>
inline Matrix<T, R, C>
operator-(const Matrix<T, R, C>& lhs, const Matrix<T, R, C>& rhs) noexcept
{
    Matrix<T, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i) {retval[i] = lhs[i] - rhs[i];}
    return retval;
}

template<typename T, typename U, std::size_t R, std::size_t C>
inline Matrix<T, R, C>
operator*(const Matrix<T, R, C>& lhs, const U& rhs) noexcept
{
    static_assert(std::is_convertible<U, T>::value, "");
    Matrix<T, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i) {retval[i] = lhs[i] * rhs;}
    return retval;
}
template<typename T, typename U, std::size_t R, std::size_t C>
inline Matrix<T, R, C>
operator*(const U& lhs, const Matrix<T, R, C>& rhs) noexcept
{
    static_assert(std::is_convertible<U, T>::value, "");
    Matrix<T, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i) {retval[i] = lhs * rhs[i];}
    return retval;
}

template<typename T, typename U, std::size_t R, std::size_t C>
inline Matrix<T, R, C>
operator/(const Matrix<T, R, C>& lhs, const U& rhs) noexcept
{
    static_assert(std::is_convertible<U, T>::value, "");
    Matrix<T, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i) {retval[i] = lhs[i] / rhs;}
    return retval;
}

template<typename T, std::size_t L, std::size_t M, std::size_t N>
inline Matrix<T, L, N>
operator*(const Matrix<T, L, M>& lhs, const Matrix<T, M, N>& rhs) noexcept
{
    Matrix<T, L, N> retval;
    for(std::size_t i=0; i < L; ++i)
    {
        for(std::size_t j=0; j < N; ++j)
        {
            for(std::size_t k=0; k < M; ++k)
            {
                retval(i, j) += lhs(i, k) * rhs(k, j);
            }
        }
    }
    return retval;
}

template<typename T, std::size_t R, std::size_t C>
Matrix<T, C, R> transpose(const Matrix<T, R, C>& mat) noexcept
{
    Matrix<T, C, R> retval;
    for(std::size_t i=0; i<R; ++i)
    {
        for(std::size_t j=0; j<C; ++j)
        {
            retval(j, i) = mat(i, j);
        }
    }
    return retval;
}

// ---------------------------------------------------------------------------
// math functions...
// ---------------------------------------------------------------------------

template<typename T, std::size_t R, std::size_t C>
inline Matrix<T, R, C>
min(const Matrix<T, R, C>& lhs, const Matrix<T, R, C>& rhs) noexcept
{
    Matrix<T, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i) {retval[i] = std::min(lhs[i], rhs[i]);}
    return retval;
}
template<typename T, std::size_t R, std::size_t C>
inline Matrix<T, R, C>
max(const Matrix<T, R, C>& lhs, const Matrix<T, R, C>& rhs) noexcept
{
    Matrix<T, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i) {retval[i] = std::max(lhs[i], rhs[i]);}
    return retval;
}

template<typename T, std::size_t R, std::size_t C>
inline Matrix<T, R, C> floor(const Matrix<T, R, C>& rhs) noexcept
{
    Matrix<T, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i) {retval[i] = std::floor(rhs[i]);}
    return retval;
}
template<typename T, std::size_t R, std::size_t C>
inline Matrix<T, R, C> ceil(const Matrix<T, R, C>& rhs) noexcept
{
    Matrix<T, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i) {retval[i] = std::ceil(rhs[i]);}
    return retval;
}

// fmadd(a, b, c) := a * b + c
template<typename T, std::size_t R, std::size_t C>
inline Matrix<T, R, C>
fmadd(const T a, const Matrix<T, R, C>& b, const Matrix<T, R, C>& c) noexcept
{
    Matrix<T, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i)
    {
        retval[i] = std::fma(a, b[i], c[i]);
    }
    return retval;
}

// fmsub(a, b, c) := a * b - c
template<typename T, std::size_t R, std::size_t C>
inline Matrix<T, R, C>
fmsub(const T a, const Matrix<T, R, C>& b, const Matrix<T, R, C>& c) noexcept
{
    Matrix<T, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i)
    {
        retval[i] = std::fma(a, b[i], -c[i]);
    }
    return retval;
}

// fnmadd(a, b, c) := -a * b + c
template<typename T, std::size_t R, std::size_t C>
inline Matrix<T, R, C>
fnmadd(const T a, const Matrix<T, R, C>& b, const Matrix<T, R, C>& c) noexcept
{
    Matrix<T, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i)
    {
        retval[i] = std::fma(-a, b[i], c[i]);
    }
    return retval;
}

// fnmsub(a, b, c) := -a * b - c
template<typename T, std::size_t R, std::size_t C>
inline Matrix<T, R, C>
fnmsub(const T a, const Matrix<T, R, C>& b, const Matrix<T, R, C>& c) noexcept
{
    Matrix<T, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i)
    {
        retval[i] = -std::fma(a, b[i], c[i]);
    }
    return retval;
}

// ---------------------------------------------------------------------------
// the following functions are only for 3*3 matrices ...
// ---------------------------------------------------------------------------

template<typename T>
inline T determinant(const Matrix<T, 3, 3>& mat) noexcept
{
    return mat(0, 0) * mat(1, 1) * mat(2, 2) +
           mat(1, 0) * mat(2, 1) * mat(0, 2) +
           mat(2, 0) * mat(0, 1) * mat(1, 2) -
           mat(0, 0) * mat(2, 1) * mat(1, 2) -
           mat(2, 0) * mat(1, 1) * mat(0, 2) -
           mat(1, 0) * mat(0, 1) * mat(2, 2);
}

template<typename realT>
Matrix<realT, 3, 3>
inverse(const Matrix<realT, 3, 3>& mat, const realT det) noexcept
{
    const auto det_inv = realT(1) / det;

    Matrix<realT, 3, 3> inv;
    inv(0,0) = det_inv * (mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1));
    inv(1,1) = det_inv * (mat(0,0) * mat(2,2) - mat(0,2) * mat(2,0));
    inv(2,2) = det_inv * (mat(0,0) * mat(1,1) - mat(0,1) * mat(1,0));

    inv(0,1) = det_inv * (mat(0,2) * mat(2,1) - mat(0,1) * mat(2,2));
    inv(0,2) = det_inv * (mat(0,1) * mat(1,2) - mat(0,2) * mat(1,1));
    inv(1,2) = det_inv * (mat(0,2) * mat(1,0) - mat(0,0) * mat(1,2));

    inv(1,0) = det_inv * (mat(1,2) * mat(2,0) - mat(1,0) * mat(2,2));
    inv(2,0) = det_inv * (mat(1,0) * mat(2,1) - mat(2,0) * mat(1,1));
    inv(2,1) = det_inv * (mat(2,0) * mat(0,1) - mat(0,0) * mat(2,1));
    return inv;
}

template<typename realT>
inline Matrix<realT, 3, 3> inverse(const Matrix<realT, 3, 3>& mat) noexcept
{
    return inverse(mat, determinant(mat));
}

} // math
} // mjolnir
#endif /* MJOLNIR_MATH_MATRIX */
