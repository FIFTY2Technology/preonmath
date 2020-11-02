/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

namespace Preon
{
namespace Math
{
    // Constructor

    template<size_t M, size_t N, typename T>
    matrix<M, N, T>::matrix(const std::initializer_list<T>& values)
    {
        THROW_EXCEPTION(values.size() != M * N, std::length_error("The length of the intializer list does not match the size of the matrix"));
        size_t i = 0;
        for (const T& elem : values)
        {
            size_t column = i % N;
            size_t row = i / N;
            element(row, column) = elem;
            ++i;
        }
    }

    template<size_t M, size_t N, typename T>
    template<typename Lambda>
    matrix<M, N, T>::matrix(const Lambda& lambda, enable_if_not_convertible<Lambda, T>)
    {
        for (size_t c = 0; c < N; c++)
            for (size_t r = 0; r < M; r++)
                m_Data[c][r] = lambda(r, c);
    }

    // Column getter & setter

    template<size_t M, size_t N, typename T>
    inline vec<M, T>& matrix<M, N, T>::column(size_t column)
    {
        return m_Data[column];
    }

    template<size_t M, size_t N, typename T>
    inline void matrix<M, N, T>::setColumn(size_t column, const vec<M, T>& value)
    {
        m_Data[column] = value;
    }

    // Row getter & setter

    template<size_t M, size_t N, typename T>
    inline vec<N, T> matrix<M, N, T>::row(size_t row) const
    {
        vec<N, T> value;
        for (size_t d = 0; d < N; d++)
            value[d] = m_Data[d][row];
        return value;
    }

    template<size_t M, size_t N, typename T>
    inline void matrix<M, N, T>::setRow(size_t row, const vec<N, T>& value)
    {
        for (size_t d = 0; d < N; d++)
            m_Data[d][row] = value[d];
    }

    // Filling

    template<size_t M, size_t N, typename T>
    void matrix<M, N, T>::fill(T value)
    {
        for (size_t c = 0; c < N; c++)
            for (size_t r = 0; r < M; r++)
                m_Data[c][r] = value;
    }

    template<size_t M, size_t N, typename T>
    template<size_t _M, size_t _N>
    void matrix<M, N, T>::setToIdentity(typename std::enable_if<_M == _N, T>::type*)
    {
        setToZero();

        for (size_t d = 0; d < M; d++)
            m_Data[d][d] = Simd::scalar_cast<float, T>(1);
    }

    // Transpose

    template<size_t M, size_t N, typename T>
    template<size_t _M, size_t _N>
    void matrix<M, N, T>::transpose(typename std::enable_if<_M == _N, T>::type*)
    {
        // Idea: Only go over the lower left corner of the matrix and swap the values.
        for (size_t r = 1; r < M; ++r)
            for (size_t c = 0; c < r; ++c)
                std::swap(element(r, c), element(c, r));
    }

    template<size_t M, size_t N, typename T>
    matrix<N, M, T> matrix<M, N, T>::transposed() const
    {
        matrix<N, M, T> mt;
        // Each column in m becomes a row in mt.
        for (size_t c = 0; c < N; ++c)
            mt.setRow(c, this->column(c));
        return mt;
    }

    // Inverse

    template<size_t M, size_t N, typename T>
    template<size_t _M, size_t _N>
    matrix<N, M, T> matrix<M, N, T>::inverted(typename std::enable_if<_M == _N && _M == 1, T>::type*) const
    {
        // THROW_EXCEPTION(IsZero(element(0, 0)), std::exception("matrix not invertible"));
        if (IsZero(element(0, 0)))
            return matrix<N, M, T>::identity();

        return matrix<N, M, T>(1 / element(0, 0));
    }

    template<size_t M, size_t N, typename T>
    template<size_t _M, size_t _N>
    matrix<N, M, T> matrix<M, N, T>::inverted(typename std::enable_if<_M == _N && _M == 2, T>::type*) const
    {
        matrix<N, M, T> inv;
        T det = this->determinant();
        // THROW_EXCEPTION(IsZero(det), std::exception("matrix not invertible"));
        if (IsZero(det))
            return matrix<N, M, T>::identity();
        T invDet = 1 / det;

        inv(0, 0) = invDet * element(1, 1);
        inv(0, 1) = -invDet * element(0, 1);
        inv(1, 0) = -invDet * element(1, 0);
        inv(1, 1) = invDet * element(0, 0);

        return inv;
    }

    template<typename T>
    T determinant_m22(const matrix<3, 3, T>& m, int row0, int col0, int row1, int col1)
    {
        // ad - bc for {[a b], [c d]} matrix.
        return m(row0, col0) * m(row1, col1) - m(row0, col1) * m(row1, col0);
    }

    template<size_t M, size_t N, typename T>
    template<size_t _M, size_t _N>
    matrix<N, M, T> matrix<M, N, T>::inverted(typename std::enable_if<_M == _N && _M == 3, T>::type*) const
    {
        T det = this->determinant();
        // THROW_EXCEPTION(IsZero(det), std::exception("matrix not invertible"));
        if (IsZero(det))
            return matrix<N, M, T>::identity();
        T invDet = 1 / det;
        return invDet * adjugate();  // Inverse is just the inverse determinant times the transposed cofactor matrix (which is again the adjugate).
    }

    template<typename T>
    T determinant_m33(const matrix<4, 4, T>& m, int col0, int col1, int col2, int row0, int row1, int row2)
    {
        return m(row0, col0) * (m(row1, col1) * m(row2, col2) - m(row2, col1) * m(row1, col2)) - m(row0, col1) * (m(row1, col0) * m(row2, col2) - m(row2, col0) * m(row1, col2))
            + m(row0, col2) * (m(row1, col0) * m(row2, col1) - m(row2, col0) * m(row1, col1));
    }

    template<size_t M, size_t N, typename T>
    template<size_t _M, size_t _N>
    matrix<N, M, T> matrix<M, N, T>::inverted(typename std::enable_if<_M == _N && _M == 4, T>::type*) const
    {
        matrix<N, M, T> inv;
        T det = this->determinant();
        // THROW_EXCEPTION(IsZero(det), std::exception("matrix not invertible"));
        if (IsZero(det))
            return matrix<N, M, T>::identity();
        T invDet = 1 / det;

        inv(0, 0) = determinant_m33(*this, 1, 2, 3, 1, 2, 3) * invDet;
        inv(1, 0) = -determinant_m33(*this, 0, 2, 3, 1, 2, 3) * invDet;
        inv(2, 0) = determinant_m33(*this, 0, 1, 3, 1, 2, 3) * invDet;
        inv(3, 0) = -determinant_m33(*this, 0, 1, 2, 1, 2, 3) * invDet;
        inv(0, 1) = -determinant_m33(*this, 1, 2, 3, 0, 2, 3) * invDet;
        inv(1, 1) = determinant_m33(*this, 0, 2, 3, 0, 2, 3) * invDet;
        inv(2, 1) = -determinant_m33(*this, 0, 1, 3, 0, 2, 3) * invDet;
        inv(3, 1) = determinant_m33(*this, 0, 1, 2, 0, 2, 3) * invDet;
        inv(0, 2) = determinant_m33(*this, 1, 2, 3, 0, 1, 3) * invDet;
        inv(1, 2) = -determinant_m33(*this, 0, 2, 3, 0, 1, 3) * invDet;
        inv(2, 2) = determinant_m33(*this, 0, 1, 3, 0, 1, 3) * invDet;
        inv(3, 2) = -determinant_m33(*this, 0, 1, 2, 0, 1, 3) * invDet;
        inv(0, 3) = -determinant_m33(*this, 1, 2, 3, 0, 1, 2) * invDet;
        inv(1, 3) = determinant_m33(*this, 0, 2, 3, 0, 1, 2) * invDet;
        inv(2, 3) = -determinant_m33(*this, 0, 1, 3, 0, 1, 2) * invDet;
        inv(3, 3) = determinant_m33(*this, 0, 1, 2, 0, 1, 2) * invDet;

        return inv;
    }

    // Adjungate.
    template<size_t M, size_t N, typename T>
    template<size_t _M, size_t _N>
    matrix<N, M, T> matrix<M, N, T>::adjugate(typename std::enable_if<_M == _N && _M == 3, T>::type*) const
    {
        matrix<M, N, T> out;
        out(0, 0) = determinant_m22(*this, 1, 1, 2, 2);
        out(0, 1) = determinant_m22(*this, 0, 2, 2, 1);
        out(0, 2) = determinant_m22(*this, 0, 1, 1, 2);
        out(1, 0) = determinant_m22(*this, 1, 2, 2, 0);
        out(1, 1) = determinant_m22(*this, 0, 0, 2, 2);
        out(1, 2) = determinant_m22(*this, 0, 2, 1, 0);
        out(2, 0) = determinant_m22(*this, 1, 0, 2, 1);
        out(2, 1) = determinant_m22(*this, 0, 1, 2, 0);
        out(2, 2) = determinant_m22(*this, 0, 0, 1, 1);
        return out;
    }

    // Reductions

    template<size_t M, size_t N, typename T>
    template<typename E>
    enable_if_non_simd<E, E> matrix<M, N, T>::norm(E p) const
    {
        T result = 0;
        for (size_t c = 0; c < N; c++)
        {
            for (size_t r = 0; r < M; r++)
            {
                result += std::pow(m_Data[c][r], p);
            }
        }
        result = std::pow(result, 1 / p);
        return result;
    }

    // TODO: optimize it!
    template<size_t M, size_t N, typename T>
    T determinant(const matrix<M, N, T>& mat, size_t n)
    {
        THROW_EXCEPTION(n == 0, std::invalid_argument("n must be >= 1"));
        T d = 0;
        matrix<M, N, T> temp;
        if (n == 1)
        {
            return mat(0, 0);
        }
        else if (n == 2)
        {
            d = mat(0, 0) * mat(1, 1) - mat(0, 1) * mat(1, 0);
            return d;
        }
        else
        {
            for (size_t p = 0; p < n; p++)
            {
                size_t h = 0;
                size_t k = 0;
                for (size_t i = 1; i < n; i++)
                {
                    for (size_t j = 0; j < n; j++)
                    {
                        if (j != p)
                        {
                            temp(h, k) = mat(i, j);
                            k++;
                            if (k + 1 == n)  // TODO: check if this is ok, before: k == n - 1, change because of size_t and cast stuff.
                            {
                                h++;
                                k = 0;
                            }
                        }
                    }
                }
                d = d + mat(0, p) * std::pow(-1, p) * determinant(temp, n - 1);
            }
            return d;
        }
    }

    template<size_t M, size_t N, typename T>
    template<size_t _M, size_t _N>
    T matrix<M, N, T>::determinant(typename std::enable_if<_M == _N && _M != 2 && _M != 3, T>::type*) const
    {
        return ::determinant(*this, M);
    }

    template<size_t M, size_t N, typename T>
    template<size_t _M, size_t _N>
    T matrix<M, N, T>::determinant(typename std::enable_if<_M == 3 && _N == 3, T>::type*) const
    {
        return m_Data[0][0] * (m_Data[1][1] * m_Data[2][2] - m_Data[1][2] * m_Data[2][1]) - m_Data[1][0] * (m_Data[0][1] * m_Data[2][2] - m_Data[0][2] * m_Data[2][1])
            + m_Data[2][0] * (m_Data[0][1] * m_Data[1][2] - m_Data[0][2] * m_Data[1][1]);
    }

    template<size_t M, size_t N, typename T>
    template<size_t _M, size_t _N>
    T matrix<M, N, T>::determinant(typename std::enable_if<_M == 2 && _N == 2, T>::type*) const
    {
        return (*this)(0, 0) * (*this)(1, 1) - (*this)(0, 1) * (*this)(1, 0);
    }

    template<size_t M, size_t N, typename T>
    template<size_t _M, size_t _N>
    T matrix<M, N, T>::trace(typename std::enable_if<_M == _N, T>::type*) const
    {
        T trace = 0;
        for (size_t d = 0; d < M; d++)
            trace += element(d, d);

        return trace;
    }

    template<size_t M, size_t N, typename T>
    template<size_t _M, size_t _N>
    vec<M, T> matrix<M, N, T>::diagonal(typename std::enable_if<_M == _N, T>::type*) const
    {
        vec<M, T> diagonal;
        for (size_t d = 0; d < M; d++)
            diagonal[d] = element(d, d);

        return diagonal;
    }

    // Operators

    template<size_t M, size_t N, typename T>
    matrix<M, N, T>& matrix<M, N, T>::operator+=(const matrix<M, N, T>& other)
    {
        for (size_t c = 0; c < N; c++)
            this->m_Data[c] += other.m_Data[c];

        return *this;
    }

    template<size_t M, size_t N, typename T>
    matrix<M, N, T>& matrix<M, N, T>::operator-=(const matrix<M, N, T>& other)
    {
        for (size_t c = 0; c < N; c++)
            this->m_Data[c] -= other.m_Data[c];

        return *this;
    }

    template<size_t M, size_t N, typename T>
    matrix<M, N, T>& matrix<M, N, T>::operator*=(const matrix<M, N, T>& other)
    {
        matrix<M, N, T> thisTransposed = this->transposed();

        for (size_t r = 0; r < M; r++)
            for (size_t c = 0; c < N; c++)
                this->element(r, c) = thisTransposed.column(r) * other.column(c);

        return *this;
    }

    template<size_t M, size_t N, typename T>
    matrix<M, N, T>& matrix<M, N, T>::operator*=(T factor)
    {
        for (size_t c = 0; c < N; c++)
            m_Data[c] *= factor;

        return *this;
    }

    template<size_t M, size_t N, size_t O, typename T>
    matrix<M, O, T> operator*(const matrix<M, N, T>& m1, const matrix<N, O, T>& m2)
    {
        return matrix<M, O, T>([&](size_t r, size_t c) { return m1.template rowVecProduct<T>(r, m2.column(c)); });
    }

    template<typename T, size_t M, size_t N, typename T_Mat, typename T_Vec>
    vec<M, T> operator*(const matrix<M, N, T_Mat>& m, const vec<N, T_Vec>& v)
    {
        return vec<M, T>([&](size_t r) { return m.template rowVecProduct<T, T_Vec>(r, v); });
    }

    // Helpers

    template<size_t M, size_t N, typename T>
    template<typename T_Out, typename T_Vec>
    T_Out matrix<M, N, T>::rowVecProduct(size_t r, const vec<N, T_Vec>& vec) const
    {
        // We manually compute the dot product of this->row(r) and vec,
        // since calling this->row(r) would copy the data due to data layout.
        T_Out prod = element(r, 0) * vec[0];
        for (size_t c = 1; c < N; ++c)
            prod += element(r, c) * vec[c];
        return prod;
    }

    // Overloaded non-member operators

    template<size_t M, size_t N, typename T>
    bool operator==(const matrix<M, N, T>& m1, const matrix<M, N, T>& m2)
    {
        for (size_t c = 0; c < N; c++)
            for (size_t r = 0; r < M; r++)
                if (m1(r, c) != m2(r, c))
                    return false;
        return true;
    }

    template<size_t M, size_t N, typename T>
    inline enable_if_non_simd<T, std::stringstream>& operator<<(std::stringstream& stream, const matrix<M, N, T>& matrix)
    {
        stream << "[";
        for (size_t r = 0; r < M; r++)
        {
            stream << "[";
            for (size_t c = 0; c < N; c++)
            {
                stream << matrix(r, c);
                if (c < N - 1)
                    stream << ", ";
            }
            stream << "]";
            if (r < M - 1)
                stream << ", ";
        }
        stream << "]";

        return stream;
    }

    template<size_t M, size_t N, typename T>
    inline enable_if_non_simd<T, std::ostream>& operator<<(std::ostream& ostream, const matrix<M, N, T>& matrix)
    {
        std::stringstream ss;
        ss << matrix;
        return ostream << ss.rdbuf();
    }

}  // namespace Math
}  // namespace Preon
