#ifndef GENNORMDIST_H
#define GENNORMDIST_H

#include <cmath>

#define normConstCDF(alpha) 1/(2*std::tgamma(1/alpha))
#define normConstPDF(sigma, alpha) (alpha/sigma)*normConstCDF(alpha)

#define genNormalPDF(x, sigma, mu, alpha) normConstPDF(sigma, alpha) * exp(-1*pow(abs(x - mu)/(sigma), alpha))
#define genNormalCDF(x, sigma, mu, alpha)  0.5 + ((x < mu) ? -1:1) * normConstCDF(alpha) * tgamma_lower(1/alpha, pow(abs(x - mu)/(sigma), alpha))

#define normalCDF(x) (erfc(-x/sqrt(2))/2)
#define normalPDF(x) (exp(-0.5 * (x) * (x)) / sqrt(2.0 * M_PI))

enum DistributionType {
    standard,
    general,
    skewNormal
};

#endif // GENNORMDIST_H