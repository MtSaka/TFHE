#pragma once

#include "type-traits.hpp"
#include "param.hpp"
#include "gen.hpp"
#include "poly.hpp"

#include <array>
#include <memory>

template <class Parameter>
struct SecretKey {
   private:
    std::array<int, Parameter::n> lvl0;
    std::array<Poly<int, Parameter::N>, Parameter::k> lvl1;
    std::array<std::array<ModInt, Parameter::N << 1>, Parameter::k> transformed_lvl1;

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
        for (std::size_t i = 0; i < Parameter::k; ++i) {
            lvl1[i].transform(transformed_lvl1[i]);
        }
    }
    int& get_lvl0(std::size_t i) { return lvl0[i]; }
    int get_lvl0(std::size_t i) const { return lvl0[i]; }
    Poly<int, Parameter::N>& get_lvl1(std::size_t i) { return lvl1[i]; }
    Poly<int, Parameter::N> get_lvl1(std::size_t i) const { return lvl1[i]; }
    int& get_lvl1(std::size_t i, std::size_t j) { return lvl1[i][j]; }
    int get_lvl1(std::size_t i, std::size_t j) const { return lvl1[i][j]; }
    const std::array<ModInt, Parameter::N << 1>& get_transformed_lvl1(std::size_t i) const { return transformed_lvl1[i]; }
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
        TLWElvl0 tlwe;
        for (std::size_t i = 0; i < n(); ++i) {
            tlwe.a(i) = uniform_torus16_gen(rng);
        }
        Torus16 e = normal_torus16_gen(alpha(), rng);
        tlwe.b() = m + e;
        for (std::size_t i = 0; i < n(); ++i) {
            tlwe.b() += tlwe.a(i) * s.get_lvl0(i);
        }
        return tlwe;
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
    TLWElvl0& operator+=(const TLWElvl0& rhs) {
        for (std::size_t i = 0; i < n() + 1; ++i) {
            (*this)[i] += rhs[i];
        }
        return *this;
    }
    TLWElvl0& operator-=(const TLWElvl0& rhs) {
        for (std::size_t i = 0; i < n() + 1; ++i) {
            (*this)[i] -= rhs[i];
        }
        return *this;
    }
    TLWElvl0& operator*=(const int& rhs) {
        for (std::size_t i = 0; i < n() + 1; ++i) {
            (*this)[i] *= rhs;
        }
        return *this;
    }
    TLWElvl0& operator+=(const Torus16& rhs) {
        (*this).b() += rhs;
        return *this;
    }
    TLWElvl0& operator-=(const Torus16& rhs) {
        (*this).b() -= rhs;
        return *this;
    }
    TLWElvl0 operator-() const {
        TLWElvl0 res;
        for (std::size_t i = 0; i < n() + 1; ++i) {
            res[i] = -(*this)[i];
        }
        return res;
    }
    friend TLWElvl0 operator+(const TLWElvl0& lhs, const TLWElvl0& rhs) {
        return TLWElvl0(lhs) += rhs;
    }
    friend TLWElvl0 operator-(const TLWElvl0& lhs, const TLWElvl0& rhs) {
        return TLWElvl0(lhs) -= rhs;
    }
    friend TLWElvl0 operator+(const TLWElvl0& lhs, const Torus16& rhs) {
        return TLWElvl0(lhs) += rhs;
    }
    friend TLWElvl0 operator-(const TLWElvl0& lhs, const Torus16& rhs) {
        return TLWElvl0(lhs) -= rhs;
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
    Torus b() const noexcept { return datab; }
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
    TLWElvl1& operator+=(const TLWElvl1& rhs) {
        for (std::size_t i = 0; i < k(); ++i) {
            for (std::size_t j = 0; j < N(); ++j) {
                a(i, j) += rhs.a(i, j);
            }
        }
        b() += rhs.b();
        return *this;
    }
    TLWElvl1& operator-=(const TLWElvl1& rhs) {
        for (std::size_t i = 0; i < k(); ++i) {
            for (std::size_t j = 0; j < N(); ++j) {
                a(i, j) -= rhs.a(i, j);
            }
        }
        b() -= rhs.b();
        return *this;
    }
    friend TLWElvl1 operator+(const TLWElvl1& lhs, const TLWElvl1& rhs) {
        return TLWElvl1(lhs) += rhs;
    }
    friend TLWElvl1 operator-(const TLWElvl1& lhs, const TLWElvl1& rhs) {
        return TLWElvl1(lhs) -= rhs;
    }
};

template <class Parameter>
struct KeySwitchKey {
   private:
    std::array<TLWElvl0<Parameter>, Parameter::kN * Parameter::t * (1 << (Parameter::basebit - 1))> data;

   public:
    TLWElvl0<Parameter>& operator()(std::size_t i, std::size_t j, std::size_t k) {
        return data[(i * Parameter::t + j) * (1 << (Parameter::basebit - 1)) + k];
    }
    TLWElvl0<Parameter> operator()(std::size_t i, std::size_t j, std::size_t k) const {
        return data[(i * Parameter::t + j) * (1 << (Parameter::basebit - 1)) + k];
    }
    /*
    template <RandGen Gen>
    KeySwitchKey(const SecretKey<Parameter>& s, Gen& rng) {
        static constexpr std::size_t bit_width = std::numeric_limits<Torus>::digits;
        for (std::size_t i = 0; i < Parameter::kN; ++i) {
            for (std::size_t j = 0; j < Parameter::t; ++j) {
                for (std::size_t k = 1; k <= (1 << (Parameter::basebit - 1)); ++k) {
                    Torus t = static_cast<Torus>(k * s.get_lvl1(i / Parameter::N, i % Parameter::N) * (1u << (bit_width - (j + 1) * Parameter::basebit)));
                    (*this)(i, j, k - 1) = TLWElvl0<Parameter>::encrypt(s, torus_to_torus16(t), rng);
                }
            }
        }
    }*/

    template <RandGen Gen>
    static std::shared_ptr<KeySwitchKey> make_ptr(const SecretKey<Parameter>& s, Gen& rng) {
        auto ks = std::make_shared<KeySwitchKey>();
        static constexpr std::size_t bit_width = std::numeric_limits<Torus>::digits;
        for (std::size_t i = 0; i < Parameter::kN; ++i) {
            for (std::size_t j = 0; j < Parameter::t; ++j) {
                for (std::size_t k = 1; k <= (1 << (Parameter::basebit - 1)); ++k) {
                    Torus t = static_cast<Torus>(k * s.get_lvl1(i / Parameter::N, i % Parameter::N) * (1u << (bit_width - (j + 1) * Parameter::basebit)));
                    (*ks)(i, j, k - 1) = TLWElvl0<Parameter>::encrypt(s, torus_to_torus16(t), rng);
                }
            }
        }
        return ks;
    }
};
template <class Parameter>
void identity_key_switch(TLWElvl0<Parameter>& res, const TLWElvl1<Parameter>& tlwe, const KeySwitchKey<Parameter>& ks) {
    std::array<std::array<int, Parameter::t>, Parameter::kN> a_bar = {};
    res = {};
    res.b() = torus_to_torus16(tlwe.b());
    static constexpr Torus round_offset = 1 << (32 - (1 + Parameter::basebit * Parameter::t));
    for (std::size_t i = 0; i < Parameter::kN; ++i) {
        for (std::size_t j = Parameter::t; j-- > 0;) {
            a_bar[i][j] += (tlwe.a(i / Parameter::N, i % Parameter::N) + round_offset) >> (32 - (j + 1) * Parameter::basebit) & ((1 << Parameter::basebit) - 1);
            if (a_bar[i][j] >= 1 << (Parameter::basebit - 1)) {
                a_bar[i][j] -= 1 << (Parameter::basebit);
                if (j > 0) a_bar[i][j - 1]++;
            }
            const std::size_t k = abs(a_bar[i][j]);
            if (a_bar[i][j] > 0)
                res -= ks(i, j, k - 1);
            else if (a_bar[i][j] < 0)
                res += ks(i, j, k - 1);
        }
    }
}