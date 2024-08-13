#pragma once

#include "type-traits.hpp"
#include "param.hpp"
#include "tlwe.cpp"
#include "poly.hpp"
#include <iostream>

template <class Parameter>
struct SecretKeyTRLWE {
   private:
    std::array<Poly<int, Parameter::N>, Parameter::k> s;

   public:
    SecretKeyTRLWE() {}
    SecretKeyTRLWE(const SecretKeyTRLWE& key) : s(key.s) {}
    template <RandGen Gen>
    SecretKeyTRLWE(Gen& gen) {
        std::binomial_distribution<int> dist;
        for (auto& e1 : s) {
            for (auto& e2 : e1) e2 = dist(gen);
        }
    }
    Poly<int, Parameter::N>& operator[](std::size_t i) { return s[i]; }
    Poly<int, Parameter::N> operator[](std::size_t i) const { return s[i]; }
};

template <class Parameter>
struct TRLWE {
   public:
    static constexpr uint32_t N() noexcept { return Parameter::N; }
    static constexpr uint32_t k() noexcept { return Parameter::k; }
    static constexpr double alpha_bk() noexcept { return Parameter::alpha_bk; }

   private:
    std::array<Poly<Torus, N()>, k() + 1> data;

   public:
    TRLWE() : data({}) {}
    Poly<Torus, N()>& operator[](std::size_t i) noexcept { return data[i]; }
    Poly<Torus, N()> operator[](std::size_t i) const noexcept { return data[i]; }
    Poly<Torus, N()>& a(std::size_t i) noexcept { return (*this)[i + 1]; }
    Poly<Torus, N()> a(std::size_t i) const noexcept { return (*this)[i + 1]; }
    Torus& a(std::size_t i, std::size_t j) noexcept { return (*this)[i + 1][j]; }
    Torus a(std::size_t i, std::size_t j) const noexcept { return (*this)[i + 1][j]; }
    Poly<Torus, N()>& b() noexcept { return (*this)[0]; }
    Poly<Torus, N()> b() const noexcept { return (*this)[0]; }
    Torus& b(std::size_t i) noexcept { return (*this)[0][i]; }
    Torus b(std::size_t i) const noexcept { return (*this)[0][i]; }

    template <RandGen Gen>
    static TRLWE encrypt(const SecretKeyTRLWE<Parameter>& s, const Poly<Torus, N()>& m, Gen& rng) {
        TRLWE trlwe;
        for (std::size_t i = 0; i < k(); ++i) {
            for (std::size_t j = 0; j < N(); ++j) {
                trlwe.a(i, j) = uniform_torus_gen(rng);
            }
        }
        for (std::size_t i = 0; i < N(); ++i) {
            trlwe.b(i) = m[i]+normal_torus_gen(alpha_bk(), rng);
        }
        for (std::size_t i = 0; i < k(); ++i) {
            trlwe.b() += trlwe.a(i) * s[i];
        }

        return trlwe;
    }
    template <RandGen Gen>
    static TRLWE encrypt(const SecretKeyTRLWE<Parameter>& s, Poly<bool, N()> m, Gen& rng) {
        static constexpr Torus mu = 1 << (std::numeric_limits<Torus>::digits - 3);
        Poly<Torus, N()> t;
        for (std::size_t i = 0; i < N(); ++i) {
            t[i] = m[i] ? mu : -mu;
        }
        return encrypt(s, t, rng);
    }
    template <RandGen Gen>
    static TRLWE encrypt_zero(const SecretKeyTRLWE<Parameter>& s, Gen& rng) {
        Poly<Torus, N()> m = {};
        return encrypt(s, m, rng);
    }
    Poly<Torus, N()> decrypt_poly_torus(const SecretKeyTRLWE<Parameter>& s) const {
        Poly<Torus, N()> m = b();
        for (std::size_t i = 0; i < k(); ++i) {
            m -= a(i) * s[i];
        }
        return m;
    }
    Poly<bool, N()> decrypt_poly_bool(const SecretKeyTRLWE<Parameter>& s) const {
        Poly<Torus, N()> t = decrypt_poly_torus(s);
        Poly<bool, N()> m;
        for (std::size_t i = 0; i < N(); ++i) {
            m[i] = static_cast<signed_type_t<Torus>>(t[i]) > 0;
        }
        return m;
    }
};

template <class Parameter>
TLWElvl1<Parameter> sample_extract_index(const TRLWE<Parameter>& trlwe, std::size_t index) {
    TLWElvl1<Parameter> tlwe;
    tlwe.b() = trlwe.b(index);
    for (std::size_t j = 0; j < Parameter::k; ++j) {
        for (std::size_t i = 0; i <= index; ++i) {
            tlwe.a(j, i) = trlwe.a(j, index - i);
        }
        for (std::size_t i = index + 1; index < Parameter::N; ++i) {
            tlwe.a(j, i) = -trlwe.a(j, Parameter::N + index - i);
        }
    }
    return tlwe;
}