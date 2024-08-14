#pragma once
#include "type-traits.hpp"
#include <cassert>
namespace modint_impl {
using T = uint64_t;

using large_t = typename double_size_uint<T>::type;
static constexpr int lg = std::numeric_limits<T>::digits;
static constexpr T mod = 9223336852482686977;
static constexpr T r = 70368744177662;
static constexpr T r2 = 71776118524344324;
static constexpr T minv = 9223336852482686975;

static constexpr int get_lg() { return lg; }
static constexpr T reduce(large_t x) {
    large_t tmp = (x + static_cast<large_t>(static_cast<T>(x) * minv) * mod) >> lg;
    return tmp >= mod ? tmp - mod : tmp;
}
static constexpr T transform(large_t x) { return reduce(x * r2); }
}  // namespace modint_impl

struct BitReversal {
   private:
    static constexpr int lg = 10;
    std::array<uint16_t, 1 << 10> reversed;
    inline constexpr uint16_t reverse(uint16_t x) {
        x = ((x & 0x5555) << 1) | ((x & 0xaaaa) >> 1);
        x = ((x & 0x3333) << 2) | ((x & 0xcccc) >> 2);
        x = ((x & 0x0f0f) << 4) | ((x & 0xf0f0) >> 4);
        return (x << 8) | (x >> 8);
    }
    inline constexpr uint16_t reverse(uint16_t x, int len) { return reverse(x) >> (16 - len); }

   public:
    constexpr BitReversal() : reversed({}) {
        for (uint16_t i = 0; i < (1 << 10); ++i) {
            reversed[i] = reverse(i, lg);
        }
    }
    constexpr uint16_t get(uint16_t i) const { return reversed[i]; }
};
constexpr BitReversal bit_reversal;

struct ModInt {
    using T = uint64_t;

   private:
    using large_t = typename double_size_uint<T>::type;
    T val;

   public:
    constexpr ModInt() : val(0) {}
    constexpr ModInt(uint32_t x) : val(modint_impl::transform(x < (static_cast<large_t>(modint_impl::mod) << modint_impl::get_lg()) ? static_cast<large_t>(x) : static_cast<large_t>(x % modint_impl::mod))) {}
    constexpr ModInt(int x) : ModInt(static_cast<uint32_t>(x < 0 ? -x : x)) {
        if (x < 0 && val) val = modint_impl::mod - val;
    }
    constexpr ModInt(uint64_t x) : val(modint_impl::transform(x < (static_cast<large_t>(modint_impl::mod) << modint_impl::get_lg()) ? static_cast<large_t>(x) : static_cast<large_t>(x % modint_impl::mod))) {}
    constexpr ModInt(int64_t x) : ModInt(static_cast<uint32_t>(x < 0 ? -x : x)) {
        if (x < 0 && val) val = modint_impl::mod - val;
    }
    constexpr T get() const { return modint_impl::reduce(val); }
    static constexpr T get_mod() { return modint_impl::mod; }
    constexpr ModInt& operator++() {
        val += modint_impl::r;
        if (val >= modint_impl::mod) val -= modint_impl::mod;
        return *this;
    }
    constexpr ModInt operator++(int) {
        ModInt res = *this;
        ++*this;
        return res;
    }
    constexpr ModInt& operator--() {
        if (val < modint_impl::r) val += modint_impl::mod;
        val -= modint_impl::r;
        return *this;
    }
    constexpr ModInt operator--(int) {
        ModInt res = *this;
        --*this;
        return res;
    }
    constexpr ModInt& operator+=(const ModInt& r) {
        val += r.val;
        if (val >= modint_impl::mod) val -= modint_impl::mod;
        return *this;
    }
    constexpr ModInt& operator-=(const ModInt& r) {
        if (val < r.val) val += modint_impl::mod;
        val -= r.val;
        return *this;
    }
    constexpr ModInt& operator*=(const ModInt& r) {
        val = modint_impl::reduce(static_cast<large_t>(val) * r.val);
        return *this;
    }
    constexpr ModInt pow(uint64_t n) const {
        ModInt res = 1, tmp = *this;
        while (n) {
            if (n & 1) res *= tmp;
            tmp *= tmp;
            n >>= 1;
        }
        return res;
    }
    constexpr ModInt inv() const { return pow(modint_impl::mod - 2); }
    constexpr ModInt& operator/=(const ModInt& r) { return *this *= r.inv(); }
    constexpr friend ModInt operator+(const ModInt& l, const ModInt& r) { return ModInt(l) += r; }
    constexpr friend ModInt operator-(const ModInt& l, const ModInt& r) { return ModInt(l) -= r; }
    constexpr friend ModInt operator*(const ModInt& l, const ModInt& r) { return ModInt(l) *= r; }
    constexpr friend ModInt operator/(const ModInt& l, const ModInt& r) { return ModInt(l) /= r; }
    constexpr friend bool operator==(const ModInt& l, const ModInt& r) { return l.val == r.val; }
    constexpr friend bool operator!=(const ModInt& l, const ModInt& r) { return l.val != r.val; }
};

struct NthRoot {
   private:
    static constexpr unsigned int lg = 10;
    std::array<ModInt, lg + 1> root, inv_root;
    static constexpr ModInt primitive_root = ModInt(5);

   public:
    constexpr NthRoot() : root{}, inv_root{} {
        root[lg] = primitive_root.pow((ModInt::get_mod() - 1) >> lg);
        inv_root[lg] = root[lg].pow(ModInt::get_mod() - 2);
        for (int i = lg - 1; i >= 0; --i) {
            root[i] = root[i + 1] * root[i + 1];
            inv_root[i] = inv_root[i + 1] * inv_root[i + 1];
        }
    }
    static constexpr unsigned int get_lg() { return lg; }
    constexpr ModInt get(int n) const { return root[n]; }
    constexpr ModInt inv(int n) const { return inv_root[n]; }
};
constexpr NthRoot nth_root;

template <std::size_t sz>
void ntt(std::array<ModInt, sz>& a) {
    static constexpr int lg = 10;

    static constexpr ModInt im = nth_root.get(2);
    for (std::size_t i = lg; i >= 1; i -= 2) {
        if (i == 1) {
            const ModInt w = nth_root.get(i);
            for (std::size_t j = 0; j < sz; j += (1u << 1)) {
                ModInt z = 1;
                for (std::size_t k = j; k < j + (1u << (i - 1)); ++k) {
                    const ModInt x = a[k], y = a[k + (1u << (i - 1))];
                    a[k] = x + y, a[k + (1u << (i - 1))] = (x - y) * z;
                    z *= w;
                }
            }
        } else {
            const ModInt w = nth_root.get(i);
            const std::size_t offset = 1 << (i - 2);
            for (std::size_t j = 0; j < sz; j += (1u << i)) {
                ModInt z = 1;
                for (std::size_t k = j; k < j + (1u << (i - 2)); ++k) {
                    const ModInt z2 = z * z, z3 = z2 * z;
                    const ModInt c0 = a[k], c1 = a[k + offset], c2 = a[k + offset * 2], c3 = a[k + offset * 3];
                    const ModInt c0c2 = c0 + c2, c0mc2 = c0 - c2, c1c3 = c1 + c3, c1mc3im = (c1 - c3) * im;
                    a[k] = c0c2 + c1c3;
                    a[k + offset] = (c0c2 - c1c3) * z2;
                    a[k + offset * 2] = (c0mc2 + c1mc3im) * z;
                    a[k + offset * 3] = (c0mc2 - c1mc3im) * z3;
                    z *= w;
                }
            }
        }
    }

    /*
    for (uint16_t i = 0; i < sz; ++i) {
        const uint16_t j = bit_reversal.get(i);
        if (i < j) std::swap(a[i], a[j]);
    }
    for (std::size_t i = 0; i < lg; ++i) {
        const ModInt w = nth_root.get(i + 1);
        for (std::size_t j = 0; j < sz; j += (1u << (i + 1))) {
            ModInt z = ModInt(1);
            for (std::size_t k = 0; k < (1u << i); ++k) {
                ModInt x = a[j + k], y = a[j + k + (1 << i)] * z;
                a[j + k] = x + y, a[j + k + (1 << i)] = x - y;
                z *= w;
            }
        }
    }*/
}
template <std::size_t sz, bool f = true>
void intt(std::array<ModInt, sz>& a) {
    static constexpr int lg = 10;

    /*
    for (uint16_t i = 0; i < sz; ++i) {
        const uint16_t j = bit_reversal.get(i);
        if (i < j) std::swap(a[i], a[j]);
    }
    for (std::size_t i = 0; i < lg; ++i) {
        const ModInt w = nth_root.inv(i + 1);
        for (std::size_t j = 0; j < sz; j += (1u << (i + 1))) {
            ModInt z = ModInt(1);
            for (std::size_t k = 0; k < (1u << i); ++k) {
                ModInt x = a[j + k], y = a[j + k + (1 << i)] * z;
                a[j + k] = x + y, a[j + k + (1 << i)] = x - y;
                z *= w;
            }
        }
    }*/
    static constexpr ModInt im = nth_root.inv(2);
    for (std::size_t i = 2 - (lg & 1); i <= lg; i += 2) {
        if (i == 1) {
            const ModInt w = nth_root.inv(i);
            for (std::size_t j = 0; j < sz; j += (1u << i)) {
                ModInt z = 1;
                for (std::size_t k = j; k < j + (1u << (i - 1)); ++k) {
                    const ModInt x = a[k], y = a[k + (1u << (i - 1))] * z;
                    a[k] = x + y, a[k + (1u << (i - 1))] = x - y;
                    z *= w;
                }
            }
        } else {
            const ModInt w = nth_root.inv(i);
            const std::size_t offset = 1 << (i - 2);
            for (std::size_t j = 0; j < sz; j += (1u << i)) {
                ModInt z = 1;
                for (std::size_t k = j; k < j + (1u << (i - 2)); ++k) {
                    const ModInt z2 = z * z, z3 = z2 * z;
                    const ModInt c0 = a[k], c1 = a[k + offset] * z2, c2 = a[k + offset * 2] * z, c3 = a[k + offset * 3] * z3;
                    const ModInt c0c1 = c0 + c1, c0mc1 = c0 - c1, c2c3 = c2 + c3, c2mc3im = (c2 - c3) * im;
                    a[k] = c0c1 + c2c3;
                    a[k + offset] = c0mc1 + c2mc3im;
                    a[k + offset * 2] = c0c1 - c2c3;
                    a[k + offset * 3] = c0mc1 - c2mc3im;
                    z *= w;
                }
            }
        }
    }
    if constexpr (f) {
        constexpr ModInt inv_sz = ModInt(1) / sz;
        for (auto& x : a) x *= inv_sz;
    }
}