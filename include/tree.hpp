#ifndef TREE_HPP
#define TREE_HPP

#include "bsmUtils.hpp"

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
    std::array<double, BINOMIAL_HOLDER_LENGTH(depth)> tree;
    double S, K, r, y, sigma, T;
};

template <size_t depth>
BinomialTree<depth>::BinomialTree(double S, double K, double r, double y, double sigma, double T)
    : S(S), K(K), r(r), y(y), sigma(sigma), T(T) {
    tree.fill(0.0);
}

template <size_t depth>
void BinomialTree<depth>::run() {
    // Logic to be implemented
}

// ============================= TrinomialTree =============================
template <size_t depth>
class TrinomialTree {
public:
    TrinomialTree(double S, double K, double r, double y, double sigma, double T);
    virtual void run();
    double value() const {return tree[0];}

protected:
    std::array<double, TRINOMIAL_HOLDER_LENGTH(depth)> tree;
    double S, K, r, y, sigma, T;
};

template <size_t depth>
TrinomialTree<depth>::TrinomialTree(double S, double K, double r, double y, double sigma, double T)
    : S(S), K(K), r(r), y(y), sigma(sigma), T(T) {
    tree.fill(0.0);
}

template <size_t depth>
void TrinomialTree<depth>::run() {
    // Logic to be implemented
}


// ============================= NaiveConvertibleTree =============================


template <size_t depth>
class NaiveConvertibleTree : public BinomialTree<depth> {
public:
    NaiveConvertibleTree(double S, double K, double r, double y, double sigma, double T,
                        double principal, double coupon, double cs, double conversionRatio, double conversionPrice, const std::vector<size_t>& conversionDates, const std::vector<size_t>& couponDates,
                         double callablePrice, double puttablePrice);
    void run() override;

private:
    std::array<double, depth> rn_rate; // Risk-neutral rate holder

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

// NaiveConvertibleTree implementation
template <size_t depth>
NaiveConvertibleTree<depth>::NaiveConvertibleTree(double S, double K, double r, double y, double sigma, double T,
    double principal, double coupon, double cs, double conversionRatio, double conversionPrice, const std::vector<size_t>& conversionDates, const std::vector<size_t>& couponDates,
    double callablePrice, double puttablePrice)
    : BinomialTree<depth>(S, K, r, y, sigma, T), principal(principal), coupon(coupon), cs(cs), conversionRatio(conversionRatio), conversionPrice(conversionPrice),
    conversionDates(conversionDates), couponDates(couponDates), callablePrice(callablePrice), puttablePrice(puttablePrice) {
    // Logic to be implemented
}

template <size_t depth>
void NaiveConvertibleTree<depth>::run() {
    double dt = this->T / (depth - 1);
    std::cout << dt << std::endl;
    double factor = std::exp(this->sigma * dt);
    std::cout << "Factor: " << factor << std::endl;
    double moneyness = this->conversionPrice / this->S;
    size_t convert_idx_max = std::floor(0.5 * (depth - 1 - (std::log(moneyness) / std::log(factor))));
    double base = this->S * this->conversionRatio;
    double p_u, p_d;
    std::tie(p_u, p_d) = bin_probs(this->r, this->y, this->sigma, dt);

    // Initialize the terminal nodes with initial rates
    for (size_t i = 0; i < convert_idx_max; ++i) {
        this->tree[i] = base * std::pow(factor, depth - 1 - 2 * i) + this->coupon;
        rn_rate[i] = this->r; // Initial risk-free rate for convertible nodes
    }
    for (size_t i = convert_idx_max; i < depth; ++i) {
        this->tree[i] = this->principal + this->coupon;
        rn_rate[i] = this->r + this->cs; // Credit spread added for non-convertible nodes
    }
   
    #ifdef DBG
    std::cout << "Time 0: ";
    for (size_t i = 0; i < depth; ++i) {
        std::cout << this->tree[i] << " ";
    }
    std::cout << std::endl;
    #endif
   
    for (size_t i = 1; i < depth; ++i) {
        for (size_t j = 0; j < depth - i; ++j) {
            // Calculate expected risk-neutral rate
            rn_rate[j] = p_u * rn_rate[j] + p_d * rn_rate[j + 1];
            
            // Option value calculation with discounting
            this->tree[j] = std::exp(-rn_rate[j] * dt) * (p_u * this->tree[j] + p_d * this->tree[j + 1]);
        }
        // std::cout << this->tree[0] << std::endl;

        if (is_in_vec(i, this->conversionDates)) {
            double S_t = this->S * std::pow(factor, depth - (i + 1));
            for (size_t j = 0; j < depth - i; ++j) {
                this->tree[j] = std::max(this->tree[j], S_t * this->conversionRatio);
                S_t /= factor;
            }
        }

        if (is_in_vec(i, this->couponDates)) {
            for (size_t j = 0; j < depth - i; ++j) {
                this->tree[j] += this->coupon;
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

#endif