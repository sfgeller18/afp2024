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

inline std::mt19937 generator(unsigned seed = 666) {
    if (seed = 666) {return std::mt19937(seed);}
    else {return std::mt19937(std::random_device{}());}
}

inline void alignGamma(double* samples, const double& mu, const size_t& vecSize) {
    std::mt19937 gen = generator(); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(0, 1); // Define the distribution for 0 and 1
    for (size_t i = 0; i < vecSize; ++i) {
        samples[i] *= dis(gen) == 0 ? -1 : 1; // Multiply by -1 or 1 randomly
        samples[i] += mu;
    }
}



// Generate a length n, d-dimensional N(0, 1) Sobol QRNG sample
class sobolSample {
    public:
        sobolSample(const size_t n, const size_t dim, const DistributionType& type = DistributionType::standard, const std::vector<double>& kwargs = {});
        ~sobolSample();
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

class sobolNormal : public sobolSample {
public:
    sobolNormal(const size_t n, const size_t dim, const double& mu, const double& sigma)
        : sobolSample(n, dim, standard, std::vector<double>{mu, sigma}) {} // Pass "normal" as type
};


// Initialize a length n, d-dimensional N(mu, sigma, moment) Sobol QRNG sample \sim e^(-z^moment)
class sobolGenNormal : public sobolSample {
        sobolGenNormal(const size_t n, const size_t dim, const double& moment, const double& mu, const double& sigma)
        : sobolSample(n, dim, general, std::vector<double>{moment, mu, sigma}) {} // Pass "gen" as type
};

class sobolSkewNormal : public sobolSample {
    sobolSkewNormal(const size_t n, const size_t dim, const double& mu, const double& sigma, const double& sigma2, const double& moment)
        : sobolSample(n, dim, skewNormal, std::vector<double>{mu, sigma, sigma2, moment}) {} // Pass "skewNormal" as type
};

#endif // SOBOLMKL_H