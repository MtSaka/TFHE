#pragma once
#include "type-traits.hpp"

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

struct ModInt {
    using T = uint64_t;

   private:
    using large_t = typename double_size_uint<T>::type;
    T val;

   public:
    static constexpr ModInt raw(int v) {
        ModInt res;
        res.val = v;
        return res;
    }
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
/*
struct ModInt {
   private:
    static constexpr uint64_t mod = (((1ull << 32) - 1) << 32) + 1;
    uint64_t val;

   public:
    static constexpr uint64_t get_mod() { return mod; }
    static constexpr ModInt raw(uint64_t v) {
        ModInt res;
        res.val = v;
        return res;
    }
    constexpr ModInt() : val(0) {}
    constexpr ModInt(uint64_t v) : val(v + static_cast<uint32_t>(-(v >= mod))) {}
    constexpr ModInt(int64_t v) : ModInt(static_cast<uint64_t>(v < 0 ? static_cast<uint64_t>(v) - (1ull << 32) + 1 : static_cast<uint64_t>(v))) {}
    constexpr ModInt(uint32_t v) : val(v) {}
    constexpr ModInt(int32_t v) : ModInt(static_cast<int64_t>(v)) {}
    constexpr ModInt& operator+=(const ModInt& rhs) {
        (*this).val += rhs.val;
        (*this).val += static_cast<uint32_t>(-((*this).val < rhs.val || (*this).val >= mod));
        return *this;
    }
    constexpr ModInt& operator-=(const ModInt& rhs) {
        (*this).val -= +static_cast<uint32_t>(-((*this).val < rhs.val));
        (*this).val -= rhs.val;
        return *this;
    }
    constexpr ModInt& operator*=(const ModInt& rhs) {
        const __uint128_t tmp = static_cast<__uint128_t>((*this).val) * rhs.val;
        const uint32_t x0 = static_cast<uint32_t>(tmp), x1 = static_cast<uint32_t>(tmp >> 32), x2 = static_cast<uint32_t>(tmp >> 64), x3 = static_cast<uint32_t>(tmp >> 96);
        uint64_t res = ((static_cast<uint64_t>(x1) + x2) << 32) + x0 - x2 - x3;
        const uint64_t x0x1 = static_cast<uint64_t>(tmp);
        res -= static_cast<uint32_t>(-((res > x0x1) && (x2 == 0)));
        res += static_cast<uint32_t>(-((res < x0x1) && (x2 != 0)));
        (*this) = ModInt(res);
        return *this;
    }
    constexpr friend ModInt operator+(const ModInt& lhs, const ModInt& rhs) {
        return ModInt(lhs) += rhs;
    }
    constexpr friend ModInt operator-(const ModInt& lhs, const ModInt& rhs) {
        return ModInt(lhs) -= rhs;
    }
    constexpr friend ModInt operator*(const ModInt& lhs, const ModInt& rhs) {
        return ModInt(lhs) *= rhs;
    }
    constexpr ModInt pow(uint64_t n) const {
        ModInt res = ModInt::raw(1), tmp = *this;
        while (n) {
            if (n & 1) res *= tmp;
            tmp *= tmp;
            n >>= 1;
        }
        return res;
    }
    constexpr uint64_t get() const { return val; }
};*/
struct NthRoot {
   private:
    static constexpr unsigned int lg = 10;
    std::array<ModInt, lg + 1> root, inv_root;
    std::array<ModInt, 1 << (lg - 1)> nth_root, inv_nth_root;
    static constexpr ModInt primitive_root = ModInt::raw(5);

   public:
    std::array<ModInt, lg - 1> rate2, inv_rate2;
    std::array<ModInt, lg - 2> rate3, inv_rate3;
    constexpr NthRoot() : root{}, inv_root{} {
        root[lg] = primitive_root.pow((ModInt::get_mod() - 1) >> lg);
        inv_root[lg] = root[lg].pow(ModInt::get_mod() - 2);
        for (int i = lg - 1; i >= 0; --i) {
            root[i] = root[i + 1] * root[i + 1];
            inv_root[i] = inv_root[i + 1] * inv_root[i + 1];
        }
        {
            ModInt z = 1;
            for (int i = 0; i < (1 << (lg - 1)); ++i) {
                nth_root[i] = z;
                z *= root[lg];
            }
            z = 1;
            const ModInt iv = ModInt(1) / (1 << (lg - 1));
            for (int i = 0; i < (1 << (lg - 1)); ++i) {
                inv_nth_root[i] = z * iv;
                z *= inv_root[lg];
            }
        }
        {
            ModInt prod = 1, iprod = 1;
            for (int i = 0; i <= lg - 2; ++i) {
                rate2[i] = root[i + 2] * prod;
                inv_rate2[i] = inv_root[i + 2] * iprod;
                prod *= inv_root[i + 2];
                iprod *= root[i + 2];
            }
        }
        {
            ModInt prod = 1, iprod = 1;
            for (int i = 0; i <= lg - 3; ++i) {
                rate3[i] = root[i + 3] * prod;
                inv_rate3[i] = inv_root[i + 3] * iprod;
                prod *= inv_root[i + 3];
                iprod *= root[i + 3];
            }
        }
    }
    static constexpr unsigned int get_lg() { return lg; }
    constexpr ModInt get(int n) const { return root[n]; }
    constexpr ModInt inv(int n) const { return inv_root[n]; }
    constexpr ModInt nth(int n) const { return nth_root[n]; }
    constexpr ModInt nth_inv(int n) const { return inv_nth_root[n]; }
};
constexpr NthRoot nth_root;

int ntt_cnt = 0, intt_cnt = 0;
template <std::size_t sz>
void ntt(std::array<ModInt, sz>& a) {
    ntt_cnt++;
    static constexpr int lg = 9;
    static constexpr ModInt im = nth_root.get(2);
    for (std::size_t i = 0; i < sz; ++i) {
        a[i] *= nth_root.nth(i);
    }
    for (int i = lg; i >= 1; i -= 2) {
        if (i == 1) {
            ModInt z = 1;
            for (std::size_t j = 0; j < sz; j += (1u << i)) {
                for (std::size_t k = j; k < j + (1u << (i - 1)); ++k) {
                    const ModInt x = a[k], y = a[k + (1u << (i - 1))] * z;
                    a[k] = x + y, a[k + (1u << (i - 1))] = (x - y);
                }
                z *= nth_root.rate2[std::countr_zero(~(unsigned int)(j >> i))];
            }
        } else {
            const std::size_t offset = 1 << (i - 2);
            ModInt z = ModInt(1);
            for (std::size_t j = 0; j < sz; j += (1u << i)) {
                for (std::size_t k = j; k < j + (1u << (i - 2)); ++k) {
                    const ModInt z2 = z * z, z3 = z2 * z;
                    const ModInt c0 = a[k], c1 = a[k + offset] * z, c2 = a[k + offset * 2] * z2, c3 = a[k + offset * 3] * z3;
                    const ModInt c0c2 = c0 + c2, c0mc2 = c0 - c2, c1c3 = c1 + c3, c1mc3im = (c1 - c3) * im;
                    a[k] = c0c2 + c1c3;
                    a[k + offset] = c0c2 - c1c3;
                    a[k + offset * 2] = c0mc2 + c1mc3im;
                    a[k + offset * 3] = c0mc2 - c1mc3im;
                }
                z *= nth_root.rate3[std::countr_zero(~(unsigned int)(j >> i))];
            }
        }
    }
}
template <std::size_t sz, bool f = true>
void intt(std::array<ModInt, sz>& a) {
    intt_cnt++;
    static constexpr int lg = 9;
    static constexpr ModInt im = nth_root.inv(2);
    for (std::size_t i = 2 - (lg & 1); i <= lg; i += 2) {
        if (i == 1) {
            ModInt z = 1;
            for (std::size_t j = 0; j < sz; j += (1u << i)) {
                for (std::size_t k = j; k < j + (1u << (i - 1)); ++k) {
                    const ModInt x = a[k], y = a[k + (1u << (i - 1))];
                    a[k] = x + y, a[k + (1u << (i - 1))] = (x - y) * z;
                }
                z *= nth_root.inv_rate2[std::countr_zero(~(unsigned int)(j >> i))];
            }
        } else {
            const std::size_t offset = 1 << (i - 2);
            ModInt z = 1;
            for (std::size_t j = 0; j < sz; j += (1u << i)) {
                for (std::size_t k = j; k < j + (1u << (i - 2)); ++k) {
                    const ModInt z2 = z * z, z3 = z2 * z;
                    const ModInt c0 = a[k], c1 = a[k + offset], c2 = a[k + offset * 2], c3 = a[k + offset * 3];
                    const ModInt c0c1 = c0 + c1, c0mc1 = c0 - c1, c2c3 = c2 + c3, c2mc3im = (c2 - c3) * im;
                    a[k] = c0c1 + c2c3;
                    a[k + offset] = (c0mc1 + c2mc3im) * z;
                    a[k + offset * 2] = (c0c1 - c2c3) * z2;
                    a[k + offset * 3] = (c0mc1 - c2mc3im) * z3;
                }
                z *= nth_root.inv_rate3[std::countr_zero(~(unsigned int)(j >> i))];
            }
        }
    }
    if constexpr (f) {
        for (std::size_t i = 0; i < sz; ++i) {
            a[i] *= nth_root.nth_inv(i);
        }
    }
}