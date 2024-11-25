#include "Tree.hpp"
#include <iostream> // Optional for debugging/logging

template <typename T>
inline bool is_in_vec(T val, const std::vector<T>& vec) {
    return std::find(vec.begin(), vec.end(), val) != vec.end();
}

// BinomialTree implementation
template <size_t depth>
BinomialTree<depth>::BinomialTree(double S, double K, double r, double y, double sigma, double T)
    : S(S), K(K), r(r), y(y), sigma(sigma), T(T) {
    tree.fill(0.0);
}

template <size_t depth>
void BinomialTree<depth>::run() {
    // Logic to be implemented
}

// TrinomialTree implementation
template <size_t depth>
TrinomialTree<depth>::TrinomialTree(double S, double K, double r, double y, double sigma, double T)
    : S(S), K(K), r(r), y(y), sigma(sigma), T(T) {
    tree.fill(0.0);
}

template <size_t depth>
void TrinomialTree<depth>::run() {
    // Logic to be implemented
}

// NaiveConvertibleTree implementation
template <size_t depth>
NaiveConvertibleTree<depth>::NaiveConvertibleTree(double S, double K, double r, double y, double sigma, double T,
    double principal, double coupon, double conversionRatio, double conversionPrice, const std::vector<size_t>& conversionDates, const std::vector<size_t>& couponDates,
    double callablePrice, double puttablePrice)
    : BinomialTree<depth>(S, K, r, y, sigma, T), principal(principal), coupon(coupon), conversionRatio(conversionRatio), conversionPrice(conversionPrice),
    conversionDates(conversionDates), couponDates(couponDates), callablePrice(callablePrice), puttablePrice(puttablePrice) {
    // Logic to be implemented
}

template <size_t depth>
void NaiveConvertibleTree<depth>::run() {
    std::array<double, depth> v_new;
    std::array<double, depth> itm_prob;
    std::array<double, depth> itm_prob_new;
    double dt = T / (depth - 1);
    double factor = std::exp(sigma * std::sqrt(dt));
    double moneyness = conversionPrice / S;
    size_t convert_idx_max = std::floor(0.5 * (depth + std::log(moneyness)));
    double base = S * conversionRatio;
    double p_u, p_d; 
    std::tie(p_u, p_d) = bin_probs(r, y, sigma, dt);
    for (size_t i = 0; i < convert_idx_max; ++i) {
        tree[i] = base * std::pow(factor, depth - 1 - 2 * i);
        itm_prob[i] = 1.0;
    }
    for (size_t i = convert_idx_max; i < depth; ++i) {
        tree[i] = principal;
        itm_prob[i] = 0.0;
    }

    for (size_t i = 1; i < depth; ++i) {
        for (size_t j = 0; j < i; ++j) {
            itm_prob_new[j] = p_u * itm_prob[j] + p_d * itm_prob[j + 1];
            v_new[j] = std::exp(-itm_prob_new[j] * dt) * (p_u * tree[j] + p_d * tree[j + 1]);
        }
        if (is_in_vec<size_t>(i, conversionDates)) {
            double S_t = S * std::pow(factor, (depth - (i + 1)));
            for (size_t j = 0; j < i; ++j) {
                tree[j] = std::max(v_new[j], S_t * conversionRatio);
                S_t = S_t / factor;
            }
        }
        if (is_in_vec<size_t>(i, couponDates)) {
            for (size_t j = 0; j < i; ++j) {
                tree[j] += coupon;
            }
        }
        std::copy(v_new.begin(), v_new.begin() + (depth - i), tree.begin());
        std::copy(itm_prob_new.begin(), itm_prob_new.begin() + (depth - i), itm_prob.begin());
    }

}

// Explicit template instantiations for linker (if needed)
template class BinomialTree<10>;
template class TrinomialTree<10>;
template class NaiveConvertibleTree<10>;
