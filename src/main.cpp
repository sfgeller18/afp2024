#define PY_ARRAY_UNIQUE_SYMBOL AFP

#include "sp.hpp"
#include "bsmUtils.hpp"
#include "generator.hpp"
#include "tree.hpp"
#include <ctime>
// #include "breedenLitzenberger.hpp"

int main(int argc, char** argv) {
    // take 1 param the size of the sample array
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <sample size>" << std::endl;
        return 1;
    }
    size_t sampleSize = std::stoul(argv[1]);
    generator s(sampleSize, 1, standard, {0, 1});
    s.plotSamples();
}