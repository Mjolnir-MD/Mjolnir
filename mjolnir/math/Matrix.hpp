#ifndef MJOLNIR_MATH_MATRIX
#define MJOLNIR_MATH_MATRIX
#include <mjolnir/util/type_traits.hpp>
#include <iostream>
#include <array>

namespace mjolnir
{

template<typename realT, std::size_t Row, std::size_t Col>
class Matrix
{
  public:
    using real_type = realT;
    using scalar_type = real_type;
    constexpr static std::size_t dim_row = Row;
    constexpr static std::size_t dim_col = Col;
    constexpr static std::size_t number_of_element = Row * Col;
    using container      = std::array<real_type, number_of_element>;
    using iterator       = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    template<typename T>
    using is_convertible_to_real = std::is_convertible<T, real_type>;

  public:
    Matrix() : values_{{}}{}
    ~Matrix() = default;

    template<typename ... T_args, typename std::enable_if<
        (sizeof...(T_args) == number_of_element) &&
        is_all<is_convertible_to_real, T_args...>::value, std::nullptr_t
        >::type = nullptr>
    Matrix(T_args ... args) noexcept : values_{{static_cast<real_type>(args)...}}{}

    template<typename T, typename std::enable_if<
        std::is_convertible<T, real_type>::value, std::nullptr_t>::type = nullptr>
    Matrix(const std::array<T, number_of_element>& rhs) noexcept
    {
        for(std::size_t i=0; i<number_of_element; ++i)
            (*this)[i] = static_cast<real_type>(rhs[i]);
    }

    Matrix(const Matrix& mat) = default;
    Matrix(Matrix&& mat)      = default;
    Matrix& operator=(const Matrix& mat) = default;
    Matrix& operator=(Matrix&& mat)      = default;

    template<typename T, typename std::enable_if<
        std::is_convertible<T, real_type>::value, std::nullptr_t>::type = nullptr>
    Matrix(const Matrix<T, dim_row, dim_col>& rhs) noexcept
    {
        for(std::size_t i=0; i<number_of_element; ++i)
            (*this)[i] = static_cast<real_type>(rhs[i]);
    }
    template<typename T, typename std::enable_if<
        std::is_convertible<T, real_type>::value, std::nullptr_t>::type = nullptr>
    Matrix(Matrix<T, dim_row, dim_col>&& rhs) noexcept
    {
        for(std::size_t i=0; i<number_of_element; ++i)
            (*this)[i] = static_cast<real_type>(rhs[i]);
    }

    template<typename T, typename std::enable_if<
        std::is_convertible<T, real_type>::value, std::nullptr_t>::type>
    Matrix& operator=(const Matrix<T, dim_row, dim_col>& rhs) noexcept
    {
        for(std::size_t i=0; i<number_of_element; ++i)
            (*this)[i] = static_cast<real_type>(rhs[i]);
    }
    template<typename T, typename std::enable_if<
        std::is_convertible<T, real_type>::value, std::nullptr_t>::type>
    Matrix& operator=(Matrix<T, dim_row, dim_col>&& rhs) noexcept
    {
        for(std::size_t i=0; i<number_of_element; ++i)
            (*this)[i] = static_cast<real_type>(rhs[i]);
    }

    template<typename T, class = typename std::enable_if<
        std::is_convertible<T, real_type>::value>::type>
    Matrix& operator+=(const Matrix<T, dim_row, dim_col>& mat) noexcept;

    template<typename T, class = typename std::enable_if<
        std::is_convertible<T, real_type>::value>::type>
    Matrix& operator-=(const Matrix<T, dim_row, dim_col>& mat) noexcept;

    template<typename T, class = typename std::enable_if<
        std::is_convertible<T, real_type>::value>::type>
    Matrix& operator*=(const T& scl) noexcept;

    template<typename T, class = typename std::enable_if<
        std::is_convertible<T, real_type>::value>::type>
    Matrix& operator/=(const T& scl) noexcept;

    scalar_type  at(const std::size_t i, const std::size_t j) const;
    scalar_type& at(const std::size_t i, const std::size_t j);
    scalar_type  operator()(const std::size_t i, const std::size_t j) const noexcept;
    scalar_type& operator()(const std::size_t i, const std::size_t j)       noexcept;

    scalar_type  at(const std::size_t i) const {return values_.at(i);}
    scalar_type& at(const std::size_t i)       {return values_.at(i);}
    scalar_type  operator[](const std::size_t i) const noexcept {return values_[i];}
    scalar_type& operator[](const std::size_t i)       noexcept {return values_[i];}

    iterator       begin()        noexcept {return values_.begin();}
    iterator       end()          noexcept {return values_.end();}
    const_iterator begin()  const noexcept {return values_.begin();}
    const_iterator end()    const noexcept {return values_.end();}
    const_iterator cbegin() const noexcept {return values_.cbegin();}
    const_iterator cend()   const noexcept {return values_.cend();}

  private:
    container values_;
};

template<typename realT, std::size_t R, std::size_t C>
template<typename T, class>
inline Matrix<realT, R, C>&
Matrix<realT, R, C>::operator+=(const Matrix<T, dim_row, dim_col>& mat) noexcept
{
    for(std::size_t i=0; i<R*C; ++i)
        (*this)[i] += mat[i];
    return *this;
}

template<typename realT, std::size_t R, std::size_t C>
template<typename T, class>
inline Matrix<realT, R, C>&
Matrix<realT, R, C>::operator-=(const Matrix<T, dim_row, dim_col>& mat) noexcept
{
    for(std::size_t i=0; i<R*C; ++i)
        (*this)[i] -= mat[i];
    return *this;
}

template<typename realT, std::size_t R, std::size_t C>
template<typename T, class>
inline Matrix<realT, R, C>&
Matrix<realT, R, C>::operator*=(const T& s) noexcept
{
    for(std::size_t i=0; i<R*C; ++i)
        (*this)[i] *= s;
    return *this;
}

template<typename realT, std::size_t R, std::size_t C>
template<typename T, class>
inline Matrix<realT, R, C>&
Matrix<realT, R, C>::operator/=(const T& s) noexcept
{
    for(std::size_t i=0; i<R*C; ++i)
        (*this)[i] /= s;
    return *this;
}

template<typename T, typename U, std::size_t R, std::size_t C>
inline Matrix<typename std::common_type<T, U>::type, R, C>
operator+(const Matrix<T, R, C>& lhs, const Matrix<U, R, C>& rhs) noexcept
{
    Matrix<typename std::common_type<T, U>::type, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i)
        retval[i] = lhs[i] + rhs[i];
    return retval;
}

template<typename T, typename U, std::size_t R, std::size_t C>
inline Matrix<typename std::common_type<T, U>::type, R, C>
operator-(const Matrix<T, R, C>& lhs, const Matrix<U, R, C>& rhs) noexcept
{
    Matrix<typename std::common_type<T, U>::type, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i)
        retval[i] = lhs[i] - rhs[i];
    return retval;
}

template<typename T, typename U, std::size_t R, std::size_t C, class = typename
         std::enable_if<std::is_convertible<U, T>::value>::type>
inline Matrix<T, R, C>
operator*(const Matrix<T, R, C>& lhs, const U rhs) noexcept
{
    Matrix<typename std::common_type<T, U>::type, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i)
        retval[i] = lhs[i] * rhs;
    return retval;
}

template<typename T, typename U, std::size_t R, std::size_t C, class = typename
         std::enable_if<std::is_convertible<U, T>::value>::type>
inline Matrix<T, R, C>
operator*(const U lhs, const Matrix<T, R, C>& rhs) noexcept
{
    Matrix<typename std::common_type<T, U>::type, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i)
        retval[i] = lhs * rhs[i];
    return retval;
}

template<typename T, typename U, std::size_t R, std::size_t C, class = typename
         std::enable_if<std::is_convertible<U, T>::value>::type>
inline Matrix<T, R, C>
operator/(const Matrix<T, R, C>& lhs, const U rhs) noexcept
{
    Matrix<T, R, C> retval;
    for(std::size_t i=0; i<R*C; ++i)
        retval[i] = lhs[i] / rhs;
    return retval;
}

template<typename T, typename U, std::size_t L, std::size_t M, std::size_t N,
    class = typename std::enable_if<std::is_convertible<U, T>::value>::type>
inline Matrix<typename std::common_type<T, U>::type, L, N>
operator*(const Matrix<T, L, M>& lhs, const Matrix<U, M, N>& rhs) noexcept
{
    Matrix<typename std::common_type<T, U>::type, L, N> retval;
    for(std::size_t i=0; i < L; ++i)
        for(std::size_t j=0; j < N; ++j)
            for(std::size_t k=0; k < M; ++k)
                retval(i, j) += lhs(i, k) * rhs(k, j);
    return retval;
}

template<typename realT, std::size_t R, std::size_t C>
typename Matrix<realT, R, C>::scalar_type
Matrix<realT, R, C>::at(const std::size_t i, const std::size_t j) const
{
    return this->values_.at(i * C + j);
}

template<typename realT, std::size_t R, std::size_t C>
typename Matrix<realT, R, C>::scalar_type&
Matrix<realT, R, C>::at(const std::size_t i, const std::size_t j)
{
    return this->values_.at(i * C + j);
}

template<typename realT, std::size_t R, std::size_t C>
typename Matrix<realT, R, C>::scalar_type
Matrix<realT, R, C>::operator()(const std::size_t i, const std::size_t j) const noexcept
{
    return this->values_[i * C + j];
}

template<typename realT, std::size_t R, std::size_t C>
typename Matrix<realT, R, C>::scalar_type&
Matrix<realT, R, C>::operator()(const std::size_t i, const std::size_t j) noexcept
{
    return this->values_[i * C + j];
}

template<typename realT, std::size_t R, std::size_t C>
Matrix<realT, C, R> transpose(const Matrix<realT, R, C>& mat) noexcept
{
    Matrix<realT, C, R> retval;
    for(std::size_t i=0; i<R; ++i)
        for(std::size_t j=0; j<C; ++j)
            retval(j, i) = mat(i, j);
    return retval;
}

template<typename charT, typename traits,
         typename realT, std::size_t R, std::size_t C>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const Matrix<realT, R, C>& m)
{
    os << "(";
    for(std::size_t i=0; i<R*C; ++i)
        os << m[i] << " ";
    os << ")";
    return os;
}

// for 3*3 only ...
template<typename realT>
inline realT determinant(const Matrix<realT, 3, 3>& mat) noexcept
{
    return mat(0,0) * mat(1,1) * mat(2,2) +
           mat(1,0) * mat(2,1) * mat(0,2) +
           mat(2,0) * mat(0,1) * mat(1,2) -
           mat(0,0) * mat(2,1) * mat(1,2) -
           mat(2,0) * mat(1,1) * mat(0,2) -
           mat(1,0) * mat(0,1) * mat(2,2);
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

} // mjolnir
#endif /* MJOLNIR_MATH_MATRIX */
