#ifndef TREE_HPP
#define TREE_HPP

#include "bsmUtils.hpp"

#define SIMD_ALIGNMENT 64
#define SIMDLEN 8


constexpr std::pair<double, double> bin_probs(double r, double y, double sigma, double dt) {
    double u = std::exp(sigma * std::sqrt(dt)); // Up factor
    double d = 1.0 / u;                        // Down factor (inverse of up)
    double p_up = (std::exp((r - y) * dt) - d) / (u - d); // Risk-neutral upward probability
    double p_down = 1.0 - p_up;                             // Downward probability
    return {p_up, p_down};
}

constexpr size_t int_pow(size_t base, size_t exp) {
    return (exp == 0) ? 1 : base * int_pow(base, exp - 1);
}

// For the binomial case we can work on one array of length n since the last layer has n leaves (we will use a temp array to hold the t_{n - 1} values and copy them back)
constexpr size_t BINOMIAL_HOLDER_LENGTH(size_t n) {
    return n;
}

//For trinomial tree max layer size is 2n + 1
constexpr size_t TRINOMIAL_HOLDER_LENGTH(size_t n) {
    return 2 * n + 1;
}

template <typename T>
inline bool is_in_vec(T val, const std::vector<T>& vec) {
    return std::find(vec.begin(), vec.end(), val) != vec.end();
}

// ============================= BinomialTree =============================
template <size_t depth>
class BinomialTree {
public:
    BinomialTree(double S, double K, double r, double y, double sigma, double T);
    virtual void run();
    double value() const {return tree[0];}
protected:
    alignas(SIMD_ALIGNMENT) std::array<double, BINOMIAL_HOLDER_LENGTH(depth)> tree;
    double dt, factor, p_u, p_d;
    double S, K, r, y, sigma, T;
};

template <size_t depth>
BinomialTree<depth>::BinomialTree(double S, double K, double r, double y, double sigma, double T)
    : S(S), K(K), r(r), y(y), sigma(sigma), T(T) {
    tree.fill(0.0);
    dt = T / (depth - 1);
    std::tie(p_u, p_d) = bin_probs(r, y, sigma, dt);
    factor = std::exp(sigma * dt);
}

template <size_t depth>
void BinomialTree<depth>::run() {
    // Logic to be implemented
}

// ============================= NaiveConvertibleTree =============================


template <size_t depth>
class NaiveConvertibleTree : public BinomialTree<depth> {
public:
    NaiveConvertibleTree(double S, double K, double r, double y, double sigma, double T,
                        double principal, double coupon, double cs, double conversionRatio, double conversionPrice, const std::vector<size_t>& conversionDates, const std::vector<size_t>& couponDates,
                         double callablePrice = 0.0, double puttablePrice = 0.0);
    void run() override;

private:
    void initConvertibleTree(double* rn_rate_data, double* tree_data, size_t convert_idx_max, double base);
    void iterateConvertibleTree(double* rn_rate_data, double* tree_data);

private:
    alignas(SIMD_ALIGNMENT) std::array<double, depth> rn_rate; // Risk-neutral rate holder
    bool isCallable;
    bool isPuttable;

private:
    double principal;
    double coupon;
    double cs;
    double conversionRatio;
    double conversionPrice;
    std::vector<size_t> conversionDates; // Changed to vector for safer memory handling
    std::vector<size_t> couponDates;
    double callablePrice;
    double puttablePrice;
};

template <size_t depth>
NaiveConvertibleTree<depth>::NaiveConvertibleTree(double S, double K, double r, double y, double sigma, double T,
    double principal, double coupon, double cs, double conversionRatio, double conversionPrice, const std::vector<size_t>& conversionDates, const std::vector<size_t>& couponDates,
    double callablePrice, double puttablePrice)
    : BinomialTree<depth>(S, K, r, y, sigma, T), principal(principal), coupon(coupon), cs(cs), conversionRatio(conversionRatio), conversionPrice(conversionPrice),
    conversionDates(conversionDates), couponDates(couponDates), callablePrice(callablePrice), puttablePrice(puttablePrice) {
        isCallable = callablePrice > 0.0;
        isPuttable = puttablePrice > 0.0;
    }


template <size_t depth>
inline void NaiveConvertibleTree<depth>::initConvertibleTree(double* rn_rate_data, double* tree_data, size_t convert_idx_max, double base) {
    #pragma omp simd aligned(tree_data: SIMD_ALIGNMENT) aligned(rn_rate_data: SIMD_ALIGNMENT) simdlen(SIMDLEN)
    for (size_t i = 0; i < convert_idx_max; ++i) {
        tree_data[i] = base * std::pow(this->factor, depth - 1 - 2 * i) + this->coupon;
        rn_rate_data[i] = this->r; // Initial risk-free rate for convertible nodes
    }

    #pragma omp simd aligned(tree_data: SIMD_ALIGNMENT) aligned(rn_rate_data: SIMD_ALIGNMENT) simdlen(SIMDLEN)
    for (size_t i = convert_idx_max; i < depth; ++i) {
        tree_data[i] = this->principal + this->coupon;
        rn_rate_data[i] = this->r + this->cs; // Credit spread added for non-convertible nodes
    }

    #ifdef DBG
    std::cout << "Time 0: ";
    for (size_t i = 0; i < depth; ++i) {
        std::cout << this->tree[i] << " ";
    }
    std::cout << std::endl;
    for (size_t i = 0; i < depth; ++i) {
        std::cout << this->rn_rate[i] << " ";
    }
    #endif
}

template <size_t depth>
inline void NaiveConvertibleTree<depth>::iterateConvertibleTree(double* rn_rate_data, double* tree_data) {
    for (size_t i = 1; i < depth; ++i) {
        #pragma omp simd aligned(rn_rate_data:SIMD_ALIGNMENT), aligned(tree_data:SIMD_ALIGNMENT) simdlen(SIMDLEN)
        for (size_t j = 0; j < depth - i; ++j) {
            rn_rate_data[j] = this->p_u * rn_rate_data[j] + this->p_d * rn_rate_data[j + 1];
            tree_data[j] = std::exp(-rn_rate_data[j] * this->dt) * (this->p_u * tree_data[j] + this->p_d * tree_data[j + 1]);
        }

        // Apply callable and puttable price conditions
        if (isCallable) {
            #pragma omp simd aligned(tree_data: SIMD_ALIGNMENT) simdlen(SIMDLEN)
            for (size_t j = 0; j < depth - i; ++j) {
                tree_data[j] = std::min(tree_data[j], callablePrice);
            }
        }

        if (isPuttable) {
            #pragma omp simd aligned(tree_data: SIMD_ALIGNMENT) simdlen(SIMDLEN)
            for (size_t j = 0; j < depth - i; ++j) {
                tree_data[j] = std::max(tree_data[j], puttablePrice);
            }
        }

        // Handle conversion dates (apply conversion ratio)
        if (is_in_vec(i, this->conversionDates)) {
            double S_t = this->S * std::pow(this->factor, depth - (i + 1));
            for (size_t j = 0; j < depth - i; ++j) {
                tree_data[j] = std::max(tree_data[j], S_t * this->conversionRatio);
                S_t /= this->factor;
            }
        }

        // Handle coupon dates (apply coupon to tree values)
        if (is_in_vec(i, this->couponDates)) {
            #pragma omp simd aligned(tree_data: SIMD_ALIGNMENT) simdlen(SIMDLEN)
            for (size_t j = 0; j < depth - i; ++j) {
                tree_data[j] += this->coupon;
            }
        }
        #ifdef DBG
        std::cout << "Time " << i << ": " << std::endl;
        std::cout << "Rate: ";
        for (size_t j = 0; j < depth - i; ++j) {
            std::cout << rn_rate[j] << " ";
        }
        std::cout << "\n" <<  "Value: ";
        for (size_t j = 0; j < depth - i; ++j) {
            std::cout << this->tree[j] << " ";
        }
        std::cout << std::endl;
        #endif
    }

}

template <size_t depth>
void NaiveConvertibleTree<depth>::run() {
    double moneyness = this->conversionPrice / this->S;
    size_t convert_idx_max = std::floor(0.5 * (depth - 1 - (std::log(moneyness) / std::log(this->factor))));
    double base = this->S * this->conversionRatio;

    double* rn_rate_data = rn_rate.data();
    double* tree_data = this->tree.data();
    
    this->initConvertibleTree(rn_rate_data, tree_data, convert_idx_max, base);
    this->iterateConvertibleTree(rn_rate_data, tree_data);
}


#endif