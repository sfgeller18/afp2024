#ifndef SP_HPP
#define SP_HPP

#include <cmath>
#include <random>
#include <array>
#include "generator.hpp"

template <typename T, size_t N, size_t M>
class SP : generator {
    public:
    SP() : generator() {}
    SP(size_t seed) : generator(seed) {}

    void genSamples() {
        
    }

    virtual ~SP() {}
    virtual void run() = 0;

    T operator[](size_t i) {
        return x[i];
    }

    protected:
    DistributionType type;
    std::array<std::array<T, M>, N> x;

};



#endif // SP_HPP