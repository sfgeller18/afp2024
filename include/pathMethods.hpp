#ifndef PATH_METHODS_HPP
#define PATH_METHODS_HPP

#include "generator.hpp"


class SP : generator {
    public:
    SP() : generator() {}
    SP(size_t seed) : generator(seed) {}

    virtual ~SP() {}
    virtual void run() = 0;

    T operator[](size_t i) {
        return x[i];
    }

    protected:
    DistributionType type;
    const size_t N;
    const size_t M;
    std::unique_ptr<double[]> x;
};

class OU : SP {
    public:
    OU(size_t N, size_t M, double mu, double sigma, double theta, double dt) : N(N), M(M), mu(mu), sigma(sigma), theta(theta), dt(dt) {
        x = std::make_unique<double[]>(N * M);
    }

    private:
    double mu;
    double sigma;
    double theta;
    double dt;
}