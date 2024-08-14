#pragma once

#include "type-traits.hpp"
#include "param.hpp"
#include "gen.hpp"
#include "trlwe.cpp"

template <class Parameter>
struct TRGSW {
   public:
    static constexpr uint32_t N() noexcept { return Parameter::N; }
    static constexpr uint32_t k() noexcept { return Parameter::k; }
    static constexpr uint32_t l() noexcept { return Parameter::l; }

   private:
    std::array<TRLWE<Parameter>, (k() + 1) * l()> data;

   public:
    TRLWE<Parameter>& operator[](std::size_t i) { return data[i]; }
    const TRLWE<Parameter>& operator[](std::size_t i) const { return data[i]; }

    template <RandGen Gen>
    static TRGSW encrypt(const SecretKey<Parameter>& s, const Poly<int, N()>& mu, Gen& rng) {
        static constexpr uint32_t bit_width = std::numeric_limits<Torus>::digits;
        TRGSW trgsw = {};
        for (auto& trlwe : trgsw.data) trlwe = TRLWE<Parameter>::encrypt_zero(s, rng);
        for (std::size_t i = 0; i < k(); ++i) {
            for (std::size_t j = 0; j < l(); ++j) {
                for (std::size_t k = 0; k < N(); ++k) {
                    Torus t = static_cast<Torus>(mu[k]) * (1u << (bit_width - (j + 1) * Parameter::Bgbit));
                    trgsw[i * l() + j].a(i, k) += t;
                }
            }
        }
        for (std::size_t i = 0; i < l(); ++i) {
            for (std::size_t j = 0; j < N(); ++j) {
                Torus t = static_cast<Torus>(mu[j]) * (1u << (bit_width - (i + 1) * Parameter::Bgbit));
                trgsw[k() * l() + i].b(j) += t;
            }
        }
        return trgsw;
    }
    template <RandGen Gen>
    static TRGSW encrypt(const SecretKey<Parameter>& s, const bool& mu, Gen& rng) {
        Poly<int, N()> t = {};
        t[0] = mu ? 1 : 0;
        return encrypt(s, t, rng);
    }
    template <RandGen Gen>
    static TRGSW encrypt(const SecretKey<Parameter>& s, const int& mu, Gen& rng) {
        Poly<int, N()> t = {};
        t[0] = mu;
        return encrypt(s, t, rng);
    }
};

template <class Parameter>
void decompose(std::array<Poly<int, Parameter::N>, Parameter::l>& a_bar, const Poly<Torus, Parameter::N>& a) {
    static constexpr uint32_t bit_width = 32;  // std::numeric_limits<int>::digits;
    static constexpr int round_offset = 1 << (bit_width - Parameter::l * Parameter::Bgbit - 1);
    a_bar = {};
    for (int i = Parameter::l - 1; i >= 0; --i) {
        for (std::size_t j = 0; j < Parameter::N; ++j) {
            a_bar[i][j] += ((a[j] + round_offset) >> (bit_width - Parameter::Bgbit * (i + 1))) & (Parameter::Bg - 1);
            if (a_bar[i][j] >= Parameter::Bg / 2) {
                a_bar[i][j] -= Parameter::Bg;
                if (i > 0) a_bar[i - 1][j]++;
            }
        }
    }
}

template <class Parameter>
void external_product(TRLWE<Parameter>& res, const TRGSW<Parameter>& trgsw, const TRLWE<Parameter>& trlwe) {
    std::array<std::array<std::array<ModInt, Parameter::N << 1>, Parameter::l>, Parameter::k> dec_a;
    std::array<std::array<ModInt, Parameter::N << 1>, Parameter::l> dec_b;
    for (std::size_t i = 0; i < Parameter::k; ++i) {
        std::array<Poly<int, Parameter::N>, Parameter::l> dec_a_tmp;
        decompose<Parameter>(dec_a_tmp, trlwe.a(i));
        for (std::size_t j = 0; j < Parameter::l; ++j) {
            dec_a_tmp[j].transform(dec_a[i][j]);
        }
    }
    {
        std::array<Poly<int, Parameter::N>, Parameter::l> dec_b_tmp;
        decompose<Parameter>(dec_b_tmp, trlwe.b());
        for (std::size_t j = 0; j < Parameter::l; ++j) {
            dec_b_tmp[j].transform(dec_b[j]);
        }
    }

    for (std::size_t i = 0; i < Parameter::k; ++i) {
        res.a(i) = {};
    }
    res.b() = {};

    for (std::size_t i = 0; i < Parameter::k; ++i) {
        for (std::size_t j = 0; j < Parameter::l; ++j) {
            for (std::size_t k = 0; k < Parameter::k; ++k) {
                res.a(k) += trgsw[i * Parameter::l + j].a(k) * dec_a[i][j];
            }
            res.b() += trgsw[i * Parameter::l + j].b() * dec_a[i][j];
        }
    }
    for (std::size_t i = 0; i < Parameter::l; ++i) {
        for (std::size_t j = 0; j < Parameter::k; ++j) {
            res.a(j) += trgsw[Parameter::k * Parameter::l + i].a(j) * dec_b[i];
        }
        res.b() += trgsw[Parameter::k * Parameter::l + i].b() * dec_b[i];
    }

}

template <class Parameter>
void cmux(TRLWE<Parameter>& res, const TRGSW<Parameter>& cond, const TRLWE<Parameter>& thn, const TRLWE<Parameter>& els) {
    const TRLWE<Parameter>&trlwe0 = els, trlwe1 = thn;
    external_product(res, cond, trlwe1 - trlwe0);
    res += trlwe0;
}