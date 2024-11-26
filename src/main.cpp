#define PY_ARRAY_UNIQUE_SYMBOL AFP

#include "sp.hpp"
#include "bsmUtils.hpp"
#include "generator.hpp"
#include "tree.hpp"
#include <ctime>
// #include "breedenLitzenberger.hpp"

constexpr double S = 5000;
constexpr double K = 5000;
constexpr double r = 0.049;
constexpr double y = 0.011;
constexpr double sigma = 0.155;
constexpr double T = 5;
constexpr double principal = 1000;
constexpr double coupon = 50;
constexpr double cs = 0.01;
constexpr double conversionRatio = 0.2;

constexpr size_t N = 20000;

int main(int argc, char** argv) {
    // take 1 param the size of the sample array
    size_t y_length = size_t(N / T);
    std::vector<size_t> conversionDates = {y_length, 2 * y_length, 3 * y_length, 4 * y_length, 5 * y_length};
    std::vector<size_t> couponDates = {y_length, 2 * y_length, 3 * y_length, 4 * y_length, 5 * y_length};

    NaiveConvertibleTree<N> cdtree(S, K, r, y, sigma, T, principal, coupon, cs, conversionRatio, 5000, conversionDates, couponDates);
    cdtree.run();
    std::cout << cdtree.value() << std::endl;
    return 0;
}