#pragma once

#include "type-traits.hpp"
#include "ntt.hpp"

template <typename T, std::size_t sz>
struct Poly : std::array<T, sz> {
    using std::array<T, sz>::array;

    Poly& operator+=(const Poly<T, sz>& rhs) {
        for (std::size_t i = 0; i < sz; ++i) {
            (*this)[i] += rhs[i];
        }
        return *this;
    }
    Poly& operator-=(const Poly<T, sz>& rhs) {
        for (std::size_t i = 0; i < sz; ++i) {
            (*this)[i] -= rhs[i];
        }
        return *this;
    }
    Poly& operator*=(const Poly<int, sz>& rhs) {
        std::array<ModInt, sz << 1> a = {}, b = {};

        for (std::size_t i = 0; i < sz; ++i) {
            a[i] = ModInt((*this)[i]), b[i] = ModInt(rhs[i]);
        }
        // for (std::size_t i = sz; i < (sz << 1); ++i) a[i] = ModInt(), b[i] = ModInt();
        ntt(a), ntt(b);
        for (std::size_t i = 0; i < (sz << 1); ++i) {
            a[i] *= b[i];
        }
        intt(a);
        for (std::size_t i = 0; i < sz; ++i) {
            (*this)[i] = static_cast<T>(static_cast<T>(a[i].get()) - static_cast<T>(a[i + sz].get()));
        }
        return *this;
    }
    friend Poly<T, sz> operator+(const Poly<T, sz>& lhs, const Poly<T, sz>& rhs) {
        return Poly(lhs) += rhs;
    }
    friend Poly<T, sz> operator-(const Poly<T, sz>& lhs, const Poly<T, sz>& rhs) {
        return Poly(lhs) -= rhs;
    }
    friend Poly<T, sz> operator*(const Poly<T, sz>& lhs, const Poly<int, sz>& rhs) {
        return Poly(lhs) *= rhs;
    }
    Poly<T, sz> mult_pow_x_k(std::size_t k) const {
        Poly<T, sz> res;
        if (k < sz) {
            for (std::size_t i = 0; i < sz - k; ++i) {
                res[i + k] = (*this)[i];
            }
            for (std::size_t i = sz - k; i < sz; ++i) {
                res[i + k - sz] = static_cast<T>(-(*this)[i]);
            }
        } else {
            const std::size_t l = k - sz;
            for (std::size_t i = 0; i < sz - l; ++i) {
                res[i + l] = static_cast<T>(-(*this)[i]);
            }
            for (std::size_t i = sz - l; i < sz; ++i) {
                res[i + l - sz] = (*this)[i];
            }
        }
        return res;
    }
    /*
    std::array<ModInt, sz << 1> ntt() const {
        std::array<ModInt, sz << 1> a = {};
        for (std::size_t i = 0; i < sz; ++i) a[i] = (*this)[i];
        ntt(a);
        return a;
    }*/
};
