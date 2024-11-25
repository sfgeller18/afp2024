#ifndef SP_HPP
#define SP_HPP

#include <cmath>
#include <random>
#include <array>
#include "sobol.hpp"

template <typename T, size_t N>
class SP {
    public:
    SP() : gen(rd()), dis(0, 1) {}
    SP(size_t seed) : gen(seed), dis(0, 1) {}

    void genSamples() {
        for (size_t i = 0; i < N; ++i) {
            x[i] = dis(gen);
        }
    }

    virtual ~SP() {}
    virtual void run() = 0;

    T operator[](size_t i) {
        return x[i];
    }

    protected:
    std::random_device rd;
    std::mt19937 gen;
    std::normal_distribution<T> dis;
    std::array<T, N> x;

};

template <typename T, size_t N>
class OU : public SP<T, N> {
    public:
    OU() : SP<T, N>() {}
    OU(size_t seed) : SP<T, N>(seed) {}

    void run() {
        SP<T, N>::x[0] = 0.0;
        for (size_t i = 0; i < N; ++i) {
            SP<T, N>::x[i + 1] = SP<T, N>::x[i] + theta * (mu - SP<T, N>::x[i]) + sigma * SP<T, N>::dis(SP<T, N>::gen);
        }
    }

    private:
    T theta = 0.15;
    T mu = 0.0;
    T sigma = 0.2;
};

template <typename T, size_t N>
class GBM : public SP<T, N> {
    public:
    GBM() : SP<T, N>() {}
    GBM(size_t seed) : SP<T, N>(seed) {}

    void run() {
        SP<T, N>::x[0] = 1.0;
        T dt = t / N;
        for (size_t i = 0; i < N; ++i) {
            SP<T, N>::x[i + 1] = SP<T, N>::x[i] * std::exp((mu - 0.5 * sigma * sigma) * dt + sigma * std::sqrt(dt) * SP<T, N>::dis(SP<T, N>::gen));
        }
    }

    private:
    T mu = 0.1;
    T sigma = 0.2;
    T t = 1.0;
};


#endif // SP_HPP