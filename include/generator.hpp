#ifndef SOBOLMKL_H
#define SOBOLMKL_H

#include <iostream>
#include <math.h>
#include <bits/stdc++.h>
#include <string>
#include <H5Cpp.h>
// #include <numeric>
#include "genNormDist.hpp"

#include <ctime>
#include "py_manager.hpp"

#if defined(USE_MKL)
#include <mkl_vsl.h>
#include <mkl_cblas.h>
#elif defined(USE_CUDA)
#include <cuda_runtime.h>
#include <curand_kernel.h>
#else
inline std::mt19937 mtgen(unsigned seed = std::random_device{}()) {
    return std::mt19937(seed);
}

inline void alignGamma(double* samples, const double& mu, const size_t& vecSize) {
    std::mt19937 gen = mtgen(); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(0, 1); // Define the distribution for 0 and 1
    for (size_t i = 0; i < vecSize; ++i) {
        samples[i] *= dis(gen) == 0 ? -1 : 1; // Multiply by -1 or 1 randomly
        samples[i] += mu;
    }
}
#endif



// Generate a length n, d-dimensional N(0, 1) Sobol QRNG sample
class generator {
    public:
        generator(const size_t n, const size_t dim, const DistributionType& type = DistributionType::standard, const std::vector<double>& kwargs = {});
        ~generator();
        double* getSamples() const;
        double getSamples(size_t i, size_t j) const {
            return samples[i * dim + j];
        }
        void printSamples() const {
            for (size_t i = 0; i < size; ++i) {
                for (size_t j = 0; j < dim; ++j) {
                    std::cout << samples[i * dim + j] << " ";
                }
                std::cout << std::endl;
            }
        }
        void writeSamplesToFile(const std::string& filename = "") const;
        void printMoments() const;
        void plotSamples(const std::string& outputPath = "") const;
    private:
        const size_t size;
        const size_t dim;
        double* samples;
        std::vector<double> distParams;
        DistributionType distribution;
};

class qrngNormal : public generator {
public:
    qrngNormal(const size_t n, const size_t dim, const double& mu, const double& sigma)
        : generator(n, dim, standard, std::vector<double>{mu, sigma}) {} // Pass "normal" as type
};


// Initialize a length n, d-dimensional N(mu, sigma, moment) Sobol QRNG sample \sim e^(-z^moment)
class qrngGenNormal : public generator {
        qrngGenNormal(const size_t n, const size_t dim, const double& moment, const double& mu, const double& sigma)
        : generator(n, dim, general, std::vector<double>{moment, mu, sigma}) {} // Pass "gen" as type
};

class qrngSkewNormal : public generator {
    qrngSkewNormal(const size_t n, const size_t dim, const double& mu, const double& sigma, const double& sigma2, const double& moment)
        : generator(n, dim, skewNormal, std::vector<double>{mu, sigma, sigma2, moment}) {} // Pass "skewNormal" as type
};

#endif // SOBOLMKL_H