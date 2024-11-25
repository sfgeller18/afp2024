#ifndef TREE_HPP
#define TREE_HPP

#include "bsmUtils.hpp"

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

template <size_t depth>
class BinomialTree {
public:
    BinomialTree(double S, double K, double r, double y, double sigma, double T);
    virtual void run();

protected:
    std::array<double, BINOMIAL_HOLDER_LENGTH(depth)> tree;
    double S, K, r, y, sigma, T;
};

template <size_t depth>
class TrinomialTree {
public:
    TrinomialTree(double S, double K, double r, double y, double sigma, double T);
    virtual void run();

protected:
    std::array<double, TRINOMIAL_HOLDER_LENGTH(depth)> tree;
    double S, K, r, y, sigma, T;
};

template <size_t depth>
class NaiveConvertibleTree : public BinomialTree<depth> {
public:
    NaiveConvertibleTree(double S, double K, double r, double y, double sigma, double T,
                        double principal, double coupon, double conversionRatio, double conversionPrice, const std::vector<size_t>& conversionDates, const std::vector<size_t>& couponDates,
                         double callablePrice, double puttablePrice);
    void run() override;

private:
    double principal;
    double coupon;
    double conversionRatio;
    double conversionPrice;
    std::vector<size_t> conversionDates; // Changed to vector for safer memory handling
    std::vector<size_t> couponDates;
    double callablePrice;
    double puttablePrice;
};
#endif