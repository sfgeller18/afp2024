#include "sobol.hpp"


#define DEFAULT_DATA_PATH "../outputs/data/sobolSample_" + std::to_string(std::time(0)) + ".h5"
#define DEFAULT_PLOT_PATH "../outputs/plots/sobolSample_" + std::to_string(std::time(0)) + ".png"





sobolSample::sobolSample(const size_t n, const size_t dim, const DistributionType& type, const std::vector<double>& kwargs)
    : size(n), dim(dim), samples(new double[n * dim]), distribution(type) {
    std::mt19937 gen = generator();
    if (type == standard) {
        if (kwargs.size() < 2) throw std::invalid_argument("Expected mu and sigma in kwargs for standard distribution.");
        double mu = kwargs[0];
        double sigma = kwargs[1];
        distParams = {mu, sigma};

        // Generate standard normal samples and scale
        std::normal_distribution<double> dist(0.0, 1.0);
        for (size_t i = 0; i < size * dim; ++i) {
            samples[i] = dist(gen) * sigma + mu;
        }

    } else if (type == general) {
        if (kwargs.size() < 3) throw std::invalid_argument("Expected alpha, mu, and sigma in kwargs for generalized normal distribution.");
        double alpha = kwargs[0];
        double mu = kwargs[1];
        double sigma = kwargs[2];
        distParams = {alpha, mu, sigma};

        if (sigma < 1e-6) throw std::invalid_argument("Standard deviation must be positive");

        double scalingFactor = std::sqrt(std::tgamma(3 / alpha) / std::tgamma(1 / alpha)) / std::pow(sigma, alpha / 2);
        double b = std::pow(scalingFactor, -alpha);
        std::gamma_distribution<double> gammaDist(1 / alpha, b);

        // Generate samples and transform
        for (size_t i = 0; i < size * dim; ++i) {
            samples[i] = std::pow(gammaDist(gen), 1 / alpha);
        }
        alignGamma(samples, mu, size * dim);

    } else if (type == skewNormal) {
        if (kwargs.size() < 4) throw std::invalid_argument("Expected alpha, mu, sigma1, and sigma2 in kwargs for skew-normal distribution.");
        double alpha = kwargs[0];
        double mu = kwargs[1];
        double sigma1 = kwargs[2];
        double sigma2 = kwargs[3];
        distParams = {alpha, mu, sigma1, sigma2};

        if (sigma1 < 1e-6 || sigma2 < 1e-6) throw std::invalid_argument("Standard deviations must be positive");

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
    }

    // Shuffle samples for added randomness
    std::shuffle(samples, samples + size * dim, gen);
}

    

    void sobolSample::printMoments() const {
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


sobolSample::~sobolSample() {
    delete[] samples;
}

double* sobolSample::getSamples() const {
    return samples;
}

void sobolSample::writeSamplesToFile(const std::string& filename) const {
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


void sobolSample::plotSamples(const std::string& outputPath) const {
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