
#include "param.hpp"
#include "tlwe.cpp"
#include "trlwe.cpp"
#include "trgsw.cpp"
#include "bootstrapping.cpp"
#include <iostream>

template <typename Parameter>
void test_tlwe(unsigned int seed, const SecretKey<Parameter>& s) {
    std::default_random_engine rng{seed};
    {
        std::binomial_distribution<> dist;
        bool m = dist(rng);
        TLWElvl0<Parameter> tlwe = TLWElvl0<Parameter>::encrypt(s, m, rng);
        bool m_ = tlwe.decrypt_bool(s);
        std::cerr << m << " " << m_ << std::endl;
        assert(m == m_);
    }
}
template <typename Parameter>
void test_trlwe(unsigned int seed, const SecretKey<Parameter>& s) {
    std::default_random_engine rng{seed};
    {
        std::binomial_distribution<> dist;
        Poly<bool, Parameter::N> m;
        for (int i = 0; i < Parameter::N; ++i) m[i] = dist(rng);
        TRLWE<Parameter> trlwe = TRLWE<Parameter>::encrypt(s, m, rng);
        Poly<bool, Parameter::N> m_ = trlwe.decrypt_poly_bool(s);
        for (std::size_t i = 0; i < Parameter::N; ++i) {
            std::cerr << m[i] << " " << m_[i] << std::endl;
            assert(m[i] == m_[i]);
        }
        std::cerr << std::endl;
    }
}
template <typename Parameter>
void test_external_product(unsigned int seed, const SecretKey<Parameter>& s) {
    std::default_random_engine rng{seed};
    {
        std::binomial_distribution<> dist;
        Poly<bool, Parameter::N> m;
        for (int i = 0; i < Parameter::N; ++i) m[i] = dist(rng);
        int c = dist(rng) ? 1 : -1;
        auto trlwe = TRLWE<Parameter>::encrypt(s, m, rng);
        auto trgsw = TRGSW<Parameter>::encrypt(s, c, rng);
        TRLWE<Parameter> res = external_product(trgsw, trlwe);
        auto m_ = res.decrypt_poly_bool(s);
        std::cerr << c << std::endl;
        if (c == 1) {
            for (std::size_t i = 0; i < Parameter::N; ++i) {
                std::cerr << m[i] << " " << m_[i] << std::endl;
                assert(m[i] == m_[i]);
            }
        } else {
            for (std::size_t i = 0; i < Parameter::N; ++i) {
                std::cerr << m[i] << " " << m_[i] << std::endl;
                assert(m[i] != m_[i]);
            }
        }
    }
}

template <typename Parameter>
void test_cmux(unsigned int seed, const SecretKey<Parameter>& s) {
    std::default_random_engine rng{seed};
    {
        std::binomial_distribution<> dist;
        Poly<bool, Parameter::N> t = {};
        for (std::size_t i = 0; i < Parameter::N; ++i) t[i] = dist(rng);
        Poly<bool, Parameter::N> f = {};
        for (std::size_t i = 0; i < Parameter::N; ++i) f[i] = dist(rng);
        bool c = dist(rng);
        auto trlwe_t = TRLWE<Parameter>::encrypt(s, t, rng);
        auto trlwe_f = TRLWE<Parameter>::encrypt(s, f, rng);
        auto trgsw_c = TRGSW<Parameter>::encrypt(s, c, rng);
        auto trlwe = cmux(trgsw_c, trlwe_t, trlwe_f);
        Poly<bool, Parameter::N> res = trlwe.decrypt_poly_bool(s);
        for (std::size_t i = 0; i < Parameter::N; ++i) {
            std::cerr << res[i] << " " << (c ? t[i] : f[i]) << std::endl;
            assert(res[i] == (c ? t[i] : f[i]));
        }
    }
}

template <class Parameter>
void test_blind_rotate(unsigned int seed, const SecretKey<Parameter>& s, const BootstrappingKey<Parameter>& bk) {
    std::default_random_engine rng{seed};
    {
        std::binomial_distribution<> dist;
        bool m = dist(rng);
        //auto trlwe = test_vector<Parameter>();
        //trlwe=trlwe.shift(513);
        //trlwe = trlwe.shift(2);
        //trlwe = trlwe.shift(512);
        auto tlwe = TLWElvl0<Parameter>::encrypt(s, m, rng);
        TLWElvl1<Parameter> res_tlwe = gate_bootstrapping_tlwe_to_tlwe(tlwe, bk);
        bool m_ = res_tlwe.decrypt_bool(s);
        std::cerr << m << " " << m_ << std::endl;
        assert(m == m_);
    }
}
int main() {
    using P = param::Security128bit;
    const int n = 2, m = 2;
    /*
    for (int i = 0; i < n; ++i) {
        std::cerr << i << std::endl;
        SecretKey<P> key;
        {
            unsigned int seed = std::random_device{}();
            std::default_random_engine rng{seed};
            key = SecretKey<P>{rng};
        }
        for (int j = 0; j < m; ++j) {
            unsigned int seed = std::random_device{}();
            test_tlwe<P>(seed, key);
        }
    }
    for (int i = 0; i < n; ++i) {
        std::cerr << i << std::endl;
        SecretKey<P> key;
        {
            unsigned seed = std::random_device{}();
            std::default_random_engine rng{seed};
            key = SecretKey<P>{rng};
        }
        for (int j = 0; j < m; ++j) {
            unsigned int seed = std::random_device{}();
            test_trlwe<P>(seed, key);
        }
    }

    for (int i = 0; i < n; ++i) {
        std::cerr << i << std::endl;
        SecretKey<P> key;
        {
            unsigned seed = std::random_device{}();
            std::default_random_engine rng{seed};
            key = SecretKey<P>{rng};
        }
        for (int j = 0; j < m; ++j) {
            unsigned int seed = std::random_device{}();
            test_external_product<P>(seed, key);
        }
    }
    for (int i = 0; i < n; ++i) {
        std::cerr << i << std::endl;
        SecretKey<P> key;
        {
            unsigned seed = std::random_device{}();
            std::default_random_engine rng{seed};
            key = SecretKey<P>{rng};
        }
        for (int j = 0; j < m; ++j) {
            unsigned int seed = std::random_device{}();
            test_cmux<P>(seed, key);
        }
    }*/


    for (int i = 0; i < n; ++i) {
        std::cerr << i << std::endl;
        SecretKey<P> key;
        BootstrappingKey<P> bk;
        {
            unsigned seed = std::random_device{}();
            std::default_random_engine rng{seed};
            key = SecretKey<P>{rng};
            bk = BootstrappingKey<P>{key, rng};
        }
        for (int j = 0; j < m; ++j) {
            unsigned int seed = std::random_device{}();
            test_blind_rotate<P>(seed, key, bk);
        }
    }

    std::cout << "PASS" << std::endl;
}
