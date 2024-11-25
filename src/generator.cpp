#include "generator.hpp"


#define DEFAULT_DATA_PATH "../outputs/data/sobolSample_" + std::to_string(std::time(0)) + ".h5"
#define DEFAULT_PLOT_PATH "../outputs/plots/sobolSample_" + std::to_string(std::time(0)) + ".png"




// Inline function for the standard distribution
inline void generateStandardDistribution(double* samples, const size_t size, const size_t dim, const std::vector<double>& kwargs) {
#if USE_MKL
    VSLStreamStatePtr stream;
    vslNewStream(&stream, VSL_BRNG_SOBOL, static_cast<MKL_INT>(dim));
    vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, dim * size, samples, 0.0, 1.0);
    cblas_dscal(dim * size, kwargs[1], samples, 1);
    alignGamma(samples, kwargs[0], dim * size);
    vslDeleteStream(&stream);
#else
    std::mt19937 gen = mtgen();
    std::normal_distribution<double> dist(kwargs[0], kwargs[1]);
    for (size_t i = 0; i < size * dim; ++i) {
        samples[i] = dist(gen);
    }
#endif
}

// Inline function for the general distribution
inline void generateGeneralDistribution(double* samples, const size_t size, const size_t dim, const std::vector<double>& kwargs) {
#if USE_MKL
    VSLStreamStatePtr stream;
    vslNewStream(&stream, VSL_BRNG_SOBOL, static_cast<MKL_INT>(dim));
    double scalingFactor = std::sqrt(std::tgamma(3 / kwargs[0]) / std::tgamma(1 / kwargs[0])) / std::pow(kwargs[2], kwargs[0] / 2);
    double b = std::pow(scalingFactor, -kwargs[0]);
    vdRngGamma_64(VSL_RNG_METHOD_GAMMA_GNORM, stream, dim * size, samples, 1 / kwargs[0], 0.0, b);
    vdPowx_64(dim * size, samples, 1 / kwargs[0], samples);
    alignGamma(samples, kwargs[1], dim * size);
    vslDeleteStream(&stream);
#else
    std::mt19937 gen = mtgen();
    double scalingFactor = std::sqrt(std::tgamma(3 / kwargs[0]) / std::tgamma(1 / kwargs[0])) / std::pow(kwargs[2], kwargs[0] / 2);
    double b = std::pow(scalingFactor, -kwargs[0]);
    std::gamma_distribution<double> gammaDist(1 / kwargs[0], b);
    for (size_t i = 0; i < size * dim; ++i) {
        samples[i] = std::pow(gammaDist(gen), 1 / kwargs[0]);
    }
    alignGamma(samples, kwargs[1], size * dim);
#endif
}

// Inline function for the skew-normal distribution
inline void generateSkewNormalDistribution(double* samples, const size_t size, const size_t dim, const std::vector<double>& kwargs) {
#if USE_MKL
    VSLStreamStatePtr stream;
    vslNewStream(&stream, VSL_BRNG_SOBOL, static_cast<MKL_INT>(dim));
    double alpha = kwargs[0], mu = kwargs[1], sigma1 = kwargs[2], sigma2 = kwargs[3];
    double prop1 = std::pow(sigma1, alpha / 2) / (std::pow(sigma1, alpha / 2) + std::pow(sigma2, alpha / 2));
    size_t inc = static_cast<size_t>(std::floor(prop1 * size * dim));

    double scalingFactor1 = std::pow(sigma1, alpha / 2) * std::sqrt(std::tgamma(1 / alpha) / std::tgamma(3 / alpha));
    double scalingFactor2 = std::pow(sigma2, alpha / 2) * std::sqrt(std::tgamma(1 / alpha) / std::tgamma(3 / alpha));
    double muReNorm = (std::tgamma(2 / alpha) / std::tgamma(1 / alpha)) * (scalingFactor2 - scalingFactor1);

    double b1 = std::pow(scalingFactor1, alpha);
    double b2 = std::pow(scalingFactor2, alpha);

    vdRngGamma_64(VSL_RNG_METHOD_GAMMA_GNORM, stream, inc, samples, 1 / alpha, 0.0, b1);
    vdRngGamma_64(VSL_RNG_METHOD_GAMMA_GNORM, stream, dim * size - inc, samples + inc, 1 / alpha, 0.0, b2);
    vdPowx_64(dim * size, samples, 1 / alpha, samples);

    for (size_t i = 0; i < dim * size; ++i) {
        samples[i] *= (i < inc) ? -1 : 1;
        samples[i] += mu - muReNorm;
    }
    vslDeleteStream(&stream);
#else
    std::mt19937 gen = mtgen();
    double alpha = kwargs[0], mu = kwargs[1], sigma1 = kwargs[2], sigma2 = kwargs[3];
    double prop1 = std::pow(sigma1, alpha / 2) / (std::pow(sigma1, alpha / 2) + std::pow(sigma2, alpha / 2));
    size_t inc = static_cast<size_t>(std::floor(prop1 * size * dim));

    double scalingFactor1 = std::pow(sigma1, alpha / 2) * std::sqrt(std::tgamma(1 / alpha) / std::tgamma(3 / alpha));
    double scalingFactor2 = std::pow(sigma2, alpha / 2) * std::sqrt(std::tgamma(1 / alpha) / std::tgamma(3 / alpha));
    double muReNorm = (std::tgamma(2 / alpha) / std::tgamma(1 / alpha)) * (scalingFactor2 - scalingFactor1);

    double b1 = std::pow(scalingFactor1, alpha);
    double b2 = std::pow(scalingFactor2, alpha);

    std::gamma_distribution<double> gammaDist1(1 / alpha, b1);
    std::gamma_distribution<double> gammaDist2(1 / alpha, b2);

    for (size_t i = 0; i < inc; ++i) {
        samples[i] = -std::pow(gammaDist1(gen), 1 / alpha);
    }
    for (size_t i = inc; i < size * dim; ++i) {
        samples[i] = std::pow(gammaDist2(gen), 1 / alpha);
    }

    for (size_t i = 0; i < size * dim; ++i) {
        samples[i] += mu - muReNorm;
    }
#endif
}

// Generator constructor
generator::generator(const size_t n, const size_t dim, const DistributionType& type, const std::vector<double>& kwargs)
    : size(n), dim(dim), samples(new double[n * dim]), distribution(type) {
    switch (type) {
        case standard:
            generateStandardDistribution(samples, size, dim, kwargs);
            break;
        case general:
            generateGeneralDistribution(samples, size, dim, kwargs);
            break;
        case skewNormal:
            generateSkewNormalDistribution(samples, size, dim, kwargs);
            break;
        default:
            throw std::invalid_argument("Invalid distribution type.");
    }

    // Shuffle samples for added randomness
    std::mt19937 gen = mtgen();
    std::shuffle(samples, samples + size * dim, gen);
}

    

    void generator::printMoments() const {
        std::vector<double> moments(4, 0.0);

        // Calculate the mean (first moment)
        double mean = std::accumulate(samples, samples + size * dim, 0.0) / (size * dim);


        // Calculate the second moment (variance)
        std::for_each(samples, samples + size * dim, [&](double val) {
            moments[1] += (val - mean) * (val - mean);
            double diff = val - mean;
            moments[2] += diff * diff * diff;
            moments[3] += diff * diff * diff * diff;
        });

        moments[1] /= (size * dim);
        moments[2] /= (size * dim);
        moments[3] /= (size * dim);

        

        std::cout << "Mean: " << mean << std::endl;
        std::cout << "Variance: " << moments[1] << std::endl;
        std::cout << "Skewness: " << moments[2] / std::pow(moments[1], 1.5) << std::endl;
        std::cout << "Kurtosis: " << moments[3] / (moments[1] * moments[1]) - 3.0 << std::endl;
    }


generator::~generator() {
    delete[] samples;
}

double* generator::getSamples() const {
    return samples;
}

void generator::writeSamplesToFile(const std::string& filename) const {
    std::string outputFilename = filename;
    if (filename == "") {outputFilename = DEFAULT_DATA_PATH;} 
    try {
        H5::H5File file(outputFilename, H5F_ACC_TRUNC);
        hsize_t dims[2] = {size, dim}; // Rows: `size`, Columns: `dim`
        H5::DataSpace dataspace(2, dims);
        H5::DataSet dataset = file.createDataSet("samples", H5::PredType::NATIVE_DOUBLE, dataspace);
        dataset.write(samples, H5::PredType::NATIVE_DOUBLE);
        std::cout << "Data successfully written to HDF5 file: " << outputFilename << std::endl;

    } catch (H5::FileIException& e) {
        std::cerr << "HDF5 File Error: " << e.getDetailMsg() << std::endl;
    } catch (H5::DataSetIException& e) {
        std::cerr << "HDF5 Dataset Error: " << e.getDetailMsg() << std::endl;
    } catch (H5::DataSpaceIException& e) {
        std::cerr << "HDF5 Dataspace Error: " << e.getDetailMsg() << std::endl;
    } catch (...) {
        std::cerr << "An unknown error occurred while writing to HDF5 file." << std::endl;
    }
}


void generator::plotSamples(const std::string& outputPath) const {
    // Print moments
    this->printMoments();

    // Set default output filename if not provided
    std::string outputFilename = outputPath;
    if (outputPath == "") {outputFilename = DEFAULT_PLOT_PATH;}

    // Call the general Python function with NumPy array data
    npy_intp dims[] = {static_cast<npy_intp>(size * dim)};  // Flattened 1D array
    std::unordered_map<std::string, std::string> kwargs;
    kwargs["output_path"] = outputFilename;

    call_python_w_numpy(samples, dims, "plot_script", kwargs);
}