#pragma once

#include "type-traits.hpp"
#include "param.hpp"
#include "gen.hpp"
#include "poly.hpp"

#include <array>
#include <assert.h>
#include <concepts>

template <class Parameter>
struct SecretKey {
   private:
    std::array<int, Parameter::n> lvl0;
    std::array<Poly<int, Parameter::N>, Parameter::k> lvl1;

   public:
    SecretKey() {}
    SecretKey(const SecretKey<Parameter>& key) : lvl0(key.lvl0), lvl1(key.lvl1) {}
    template <RandGen Gen>
    SecretKey(Gen& rng) {
        std::binomial_distribution<int> dist;
        for (auto& e : lvl0) e = dist(rng);
        for (auto& e1 : lvl1) {
            for (auto& e2 : e1) e2 = dist(rng);
        }
    }
    int& get_lvl0(std::size_t i) { return lvl0[i]; }
    int get_lvl0(std::size_t i) const { return lvl0[i]; }
    Poly<int, Parameter::N>& get_lvl1(std::size_t i) { return lvl1[i]; }
    Poly<int, Parameter::N> get_lvl1(std::size_t i) const { return lvl1[i]; }
    int& get_lvl1(std::size_t i, std::size_t j) { return lvl1[i][j]; }
    int get_lvl1(std::size_t i, std::size_t j) const { return lvl1[i][j]; }
};

template <class Parameter>
struct TLWElvl0 {
   public:
    static constexpr std::size_t n() noexcept { return Parameter::n; }
    static constexpr double alpha() noexcept { return Parameter::alpha; }

   private:
    std::array<Torus16, n() + 1> data;

   public:
    Torus16& operator[](std::size_t i) noexcept { return data[i]; }
    Torus16 operator[](std::size_t i) const noexcept { return data[i]; }
    Torus16& a(std::size_t i) noexcept { return (*this)[i + 1]; }
    Torus16 a(std::size_t i) const noexcept { return (*this)[i + 1]; }
    Torus16& b() noexcept { return (*this)[0]; }
    Torus16 b() const noexcept { return (*this)[0]; }
    template <RandGen Gen>
    static TLWElvl0 encrypt(const SecretKey<Parameter>& s, const Torus16& m, Gen& rng) {
        TLWElvl0 TLWElvl0;
        for (std::size_t i = 0; i < n(); ++i) {
            TLWElvl0.a(i) = uniform_torus16_gen(rng);
        }
        Torus16 e = normal_torus16_gen(alpha(), rng);
        TLWElvl0.b() = m + e;
        for (std::size_t i = 0; i < n(); ++i) {
            TLWElvl0.b() += TLWElvl0.a(i) * s.get_lvl0(i);
        }
        return TLWElvl0;
    }
    template <RandGen Gen>
    static TLWElvl0 encrypt(const SecretKey<Parameter>& s, bool m, Gen& rng) {
        static constexpr Torus16 mu = 1 << (std::numeric_limits<Torus16>::digits - 3);
        return encrypt(s, static_cast<Torus16>(m ? mu : -mu), rng);
    }
    Torus16 decrypt_torus(const SecretKey<Parameter>& s) const {
        Torus16 m = b();
        for (std::size_t i = 0; i < n(); ++i) {
            m -= a(i) * s.get_lvl0(i);
        }
        return m;
    }
    bool decrypt_bool(const SecretKey<Parameter>& s) const {
        return static_cast<signed_type_t<Torus16>>(decrypt_torus(s)) > 0;
    }
};

template <class Parameter>
struct TLWElvl1 {
   public:
    static constexpr std::size_t N() noexcept { return Parameter::N; }
    static constexpr std::size_t k() noexcept { return Parameter::k; }
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
    static TLWElvl1 encrypt(const SecretKey<Parameter>& s, const Torus& m, Gen& rng) {
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
                tlwe.b() += tlwe.a(i, j) * s.get_lvl1(i, j);
            }
        }
        return tlwe;
    }
    template <RandGen Gen>
    static TLWElvl1 encrypt(const SecretKey<Parameter>& s, bool m, Gen& rng) {
        static constexpr Torus mu = 1 << (std::numeric_limits<Torus>::digits - 3);
        return encrypt(s, (m ? mu : -mu), rng);
    }
    Torus decrypt_torus(const SecretKey<Parameter>& s) const {
        Torus m = b();
        for (std::size_t i = 0; i < k(); ++i) {
            for (std::size_t j = 0; j < N(); ++j) {
                m -= a(i, j) * s.get_lvl1(i, j);
            }
        }
        return m;
    }
    bool decrypt_bool(const SecretKey<Parameter>& s) const {
        return static_cast<signed_type_t<Torus>>(decrypt_torus(s)) > 0;
    }
};