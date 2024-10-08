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
        std::array<ModInt, sz> a = {}, b = {};
        for (std::size_t i = 0; i < sz; ++i) {
            a[i] = ModInt((*this)[i]), b[i] = ModInt(rhs[i]);
        }
        ntt(a), ntt(b);
        for (std::size_t i = 0; i < sz; ++i) {
            a[i] *= b[i];
        }
        intt(a);
        for (std::size_t i = 0; i < sz; ++i) {
            const uint64_t tmp = a[i].get();
            (*this)[i] = static_cast<T>(tmp >= (1ull << 61) ? -static_cast<T>(ModInt::get_mod() - tmp) : tmp);
        }
        return *this;
    }
    Poly& operator*=(const std::array<ModInt, sz>& rhs) {
        std::array<ModInt, sz> a = {};
        for (std::size_t i = 0; i < sz; ++i) {
            a[i] = ModInt((*this)[i]);
        }
        ntt(a);
        for (std::size_t i = 0; i < sz; ++i) {
            a[i] *= rhs[i];
        }
        intt(a);
        for (std::size_t i = 0; i < sz; ++i) {
            const uint64_t tmp = a[i].get();
            (*this)[i] = static_cast<T>(tmp >= (1ull << 61) ? -static_cast<T>(ModInt::get_mod() - tmp) : tmp);
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
    friend Poly<T, sz> operator*(const Poly<T, sz>& lhs, const std::array<ModInt, sz>& rhs) {
        return Poly(lhs) *= rhs;
    }
    void mult_pow_x_k(Poly<T, sz>& res, std::size_t k) const {
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
    }
    void transform(std::array<ModInt, sz>& a) const {
        a = {};
        for (std::size_t i = 0; i < sz; ++i) a[i] = (*this)[i];
        ntt(a);
        return;
    }
    std::array<ModInt, sz> transformed() const {
        std::array<ModInt, sz> a = {};
        for (std::size_t i = 0; i < sz; ++i) a[i] = (*this)[i];
        ntt(a);
        return a;
    }
};

template <std::size_t sz>
std::array<ModInt, sz> operator*(const std::array<ModInt, sz>& lhs, const std::array<ModInt, sz>& rhs) {
    std::array<ModInt, sz> res;
    for (std::size_t i = 0; i < sz; ++i) res[i] = lhs[i] * rhs[i];
    return res;
}
template <std::size_t sz>
void add(std::array<ModInt, sz>& lhs, const std::array<ModInt, sz>& rhs) {
    for (std::size_t i = 0; i < sz; ++i) lhs[i] += rhs[i];
}