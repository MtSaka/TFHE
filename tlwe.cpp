#pragma once

#include "type-traits.hpp"
#include "param.hpp"
#include "gen.hpp"

#include <array>
#include <assert.h>
#include <concepts>

template <class Parameter>
struct SecretKeyTLWElvl0 {
   private:
    std::array<int, Parameter::n> s;

   public:
    SecretKeyTLWElvl0() {}
    SecretKeyTLWElvl0(const SecretKeyTLWElvl0<Parameter>& key) : s(key.s) {}
    template <RandGen Gen>
    SecretKeyTLWElvl0(Gen& rng) {
        std::binomial_distribution<int> dist;
        for (auto& e : s) e = dist(rng);
    }
    int& operator[](std::size_t i) { return s[i]; }
    int operator[](std::size_t i) const { return s[i]; }
};

template <class Parameter>
struct SecretKeyTLWElvl1 {
   private:
    std::array<std::array<int, Parameter::N>, Parameter::k> s;

   public:
    SecretKeyTLWElvl1() {}
    SecretKeyTLWElvl1(const SecretKeyTLWElvl1<Parameter>& key) : s(key.s) {}
    template <RandGen Gen>
    SecretKeyTLWElvl1(Gen& gen) {
        std::binomial_distribution<int> dist;
        for (auto& e1 : s) {
            for (auto& e2 : s) e2 = dist(gen);
        }
    }
    std::array<int, Parameter::N>& operator[](std::size_t i) { return s[i]; }
    std::array<int, Parameter::N> operator[](std::size_t i) const { return s[i]; }
};

template <class Parameter>
struct TLWElvl0 {
   public:
    static constexpr uint32_t n() noexcept { return Parameter::n; }
    static constexpr double alpha() noexcept { return Parameter::alpha; }

   private:
    std::array<Torus, n() + 1> data;

   public:
    Torus& operator[](std::size_t i) noexcept { return data[i]; }
    Torus operator[](std::size_t i) const noexcept { return data[i]; }
    Torus& a(std::size_t i) noexcept { return (*this)[i + 1]; }
    Torus a(std::size_t i) const noexcept { return (*this)[i + 1]; }
    Torus& b() noexcept { return (*this)[0]; }
    Torus b() const noexcept { return (*this)[0]; }
    template <RandGen Gen>
    static TLWElvl0 encrypt(const SecretKeyTLWElvl0<Parameter>& s, Torus m, Gen& rng) {
        TLWElvl0 TLWElvl0;
        for (std::size_t i = 0; i < n(); ++i) {
            TLWElvl0.a(i) = uniform_torus_gen(rng);
        }
        Torus e = normal_torus_gen(alpha(), rng);
        TLWElvl0.b() = m + e;
        for (std::size_t i = 0; i < n(); ++i) {
            TLWElvl0.b() += TLWElvl0.a(i) * s[i];
        }
        return TLWElvl0;
    }
    template <RandGen Gen>
    static TLWElvl0 encrypt(const SecretKeyTLWElvl0<Parameter>& s, bool m, Gen& rng) {
        static constexpr Torus mu = 1 << (std::numeric_limits<Torus>::digits - 3);
        return encrypt(s, (m ? mu : -mu), rng);
    }
    Torus decrypt_torus(const SecretKeyTLWElvl0<Parameter>& s) const {
        Torus m = b();
        for (std::size_t i = 0; i < n(); ++i) {
            m -= a(i) * s[i];
        }
        return m;
    }
    bool decrypt_bool(const SecretKeyTLWElvl0<Parameter>& s) const {
        return static_cast<signed_type_t<Torus>>(decrypt_torus(s)) > 0;
    }
};

template <class Parameter>
struct TLWElvl1 {
   public:
    static constexpr uint32_t N() noexcept { return Parameter::N; }
    static constexpr uint32_t k() noexcept { return Parameter::k; }
    static constexpr double alpha_bk() noexcept { return Parameter::alpha_bk; }

   private:
    std::array<std::array<Torus, N()>, k()> data;
    Torus datab;

   public:
    std::array<Torus, N()>& operator[](std::size_t i) noexcept { return data[i]; }
    std::array<Torus, N()> operator[](std::size_t i) const noexcept { return data[i]; }
    std::array<Torus, N()>& a(std::size_t i) noexcept { return (*this)[i]; }
    std::array<Torus, N()> a(std::size_t i) const noexcept { return (*this)[i]; }
    Torus& a(std::size_t i, std::size_t j) noexcept { return (*this)[i][j]; }
    Torus a(std::size_t i, std::size_t j) const noexcept { return (*this)[i][j]; }
    Torus& b() noexcept { return datab; }
    Torus b(std::size_t i) const noexcept { return datab; }
    template <RandGen Gen>
    static TLWElvl1 encrypt(const SecretKeyTLWElvl1<Parameter>& s, const Torus& m, Gen& rng) {
        TLWElvl1 tlwe;
        for (std::size_t i = 0; i < k(); ++i) {
            for (std::size_t j = 0; j < N(); ++j) {
                tlwe.a(i, j) = uniform_torus_gen(rng);
            }
        }
        Torus e = normal_torus_gen(alpha_bk(), rng);
        tlwe.b() = m + e;
        for (std::size_t i = 0; i < k(); ++i) {
            for (std::size_t j = 0; j < N(); ++j) {
                tlwe.b() += tlwe.a(i, j) * s[i][j];
            }
        }
        return tlwe;
    }
    template <RandGen Gen>
    static TLWElvl1 encrypt(const SecretKeyTLWElvl0<Parameter>& s, bool m, Gen& rng) {
        static constexpr Torus mu = 1 << (std::numeric_limits<Torus>::digits - 3);
        return encrypt(s, (m ? mu : -mu), rng);
    }
    Torus decrypt_torus(const SecretKeyTLWElvl0<Parameter>& s) const {
        Torus m = b();
        for (std::size_t i = 0; i < k(); ++i) {
            for (std::size_t j = 0; j < N(); ++j) {
                m -= a(i, j) * s[i][j];
            }
        }
        return m;
    }
    bool decrypt_bool(const SecretKeyTLWElvl0<Parameter>& s) const {
        return static_cast<signed_type_t<Torus>>(decrypt_torus(s)) > 0;
    }
};