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
    static TRGSW encrypt(const SecretKeyTRLWE<Parameter>& s, const Poly<int, N()>& mu, Gen& rng) {
        static constexpr uint32_t bit_width = std::numeric_limits<Torus>::digits;
        TRGSW trgsw = {};
        for (auto& trlwe : trgsw.data) trlwe = TRLWE<Parameter>::encrypt_zero(s, rng);
        for (std::size_t i = 0; i < k(); ++i) {
            for (std::size_t j = 0; j < l(); ++j) {
                for (std::size_t k = 0; k < N(); ++k) {
                    Torus t = static_cast<Torus>(mu[k] * (1u << (bit_width - (j + 1) * Parameter::Bgbit)));
                    trgsw[i * l() + j].a(i, k) += t;
                }
            }
        }
        for (std::size_t i = 0; i < l(); ++i) {
            for (std::size_t j = 0; j < N(); ++j) {
                Torus t = static_cast<Torus>(mu[j] * (1u << (bit_width - (i + 1) * Parameter::Bgbit)));
                trgsw[k() * l() + i].b(j) += t;
            }
        }
        /*
        for (std::size_t i = 0; i < (k() + 1) * l(); ++i) {
            for (std::size_t j = 0; j < k(); ++j) {
                std::cout << "TRGSW A " << i << " " << j << " ";
                for (std::size_t k = 0; k < N(); ++k) {
                    std::cout << trgsw[i].a(j, k) << " \n"[k == N() - 1];
                }
            }
            std::cout << "TRGSW B " << i << " " << j << " ";
            for (std::size_t j = 0;j<)
        }*/
        return trgsw;
    }
    template <RandGen Gen>
    static TRGSW encrypt(const SecretKeyTRLWE<Parameter>& s, const bool& mu, Gen& rng) {
        Poly<int, N()> t = {};
        t[0] = mu ? 1 : 0;
        return encrypt(s, t, rng);
    }
    template <RandGen Gen>
    static TRGSW encrypt(const SecretKeyTRLWE<Parameter>& s, const int& mu, Gen& rng) {
        Poly<int, N()> t = {};
        t[0] = mu;
        return encrypt(s, t, rng);
    }
};

template <class Parameter>
std::array<Poly<int, Parameter::N>, Parameter::l> decompose(const Poly<Torus, Parameter::N>& a) {
    static constexpr uint32_t bit_width = std::numeric_limits<int>::digits;
    int round_offset = 1 << (bit_width - Parameter::l * Parameter::Bgbit - 1);
    std::array<Poly<int, Parameter::N>, Parameter::l> a_bar;
    for (std::size_t i = 0; i < Parameter::l; ++i) {
        for (std::size_t j = 0; j < Parameter::N; ++j) {
            a_bar[i][j] = ((a[j] + round_offset) >> (bit_width - Parameter::Bgbit * (i + 1))) & (Parameter::Bg - 1);
        }
    }
    for (int i = Parameter::l - 1; i >= 0; --i) {
        for (std::size_t j = 0; j < Parameter::N; ++j) {
            if (a_bar[i][j] >= Parameter::Bg / 2) {
                a_bar[i][j] -= Parameter::Bg;
                if (i > 0) a_bar[i - 1][j]++;
            }
        }
    }
    return a_bar;
}

template <class Parameter>
TRLWE<Parameter> external_product(const TRGSW<Parameter>& trgsw, const TRLWE<Parameter>& trlwe) {
    std::array<std::array<Poly<int, Parameter::N>, Parameter::l>, Parameter::k> dec_a;
    std::array<Poly<int, Parameter::N>, Parameter::l> dec_b;
    for (std::size_t i = 0; i < Parameter::k; ++i) {
        dec_a[i] = decompose<Parameter>(trlwe.a(i));
    }
    dec_b = decompose<Parameter>(trlwe.b());

    TRLWE<Parameter> res;
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

    return res;
}

template <class Parameter>
TRLWE<Parameter> cmux(const TRGSW<Parameter>& cond, const TRLWE<Parameter>& thn, const TRLWE<Parameter>& els) {
    const TRLWE<Parameter>&trlwe0 = els, trlwe1 = thn;

    TRLWE<Parameter> tmp;
    for (std::size_t i = 0; i < Parameter::k; ++i) {
        tmp.a(i) = trlwe1.a(i) - trlwe0.a(i);
    }
    tmp.b() = trlwe1.b() - trlwe0.b();

    TRLWE<Parameter> res = external_product(cond, tmp);
    for (std::size_t i = 0; i < Parameter::k; ++i) {
        res.a(i) += trlwe0.a(i);
    }
    res.b() += trlwe0.b();

    return res;
}