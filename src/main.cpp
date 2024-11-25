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
constexpr double T = 1;
constexpr double principal = 1000;
constexpr double coupon = 50;
constexpr double cs = 0.01;
constexpr double conversionRatio = 0.2;

int main(int argc, char** argv) {
    // take 1 param the size of the sample array
    NaiveConvertibleTree<20000> cdtree(S, K, r, y, sigma, T, principal, coupon, cs, conversionRatio, 5000, {1, 2, 3, 4, 5}, {1, 2, 3, 4, 5}, 5000, 5000);
    cdtree.run();
    std::cout << cdtree.value() << std::endl;
    return 0;
}