#include "matrixwrapper.h"
#include <stdexcept>
#include <iostream>
#include <thread>
#include <vector>
#include <cmath> 
#include <functional>
#include "ThreadManager/ThreadManager.h"
#include <random>

double MatrixWrapper::threadScalingFactor = 1.0;

int calculateThreadsBasedOnMatrixSize(int rows, int cols, double scalingFactor, int maxThreads);


MatrixWrapper::MatrixWrapper(int rows, int cols, double value) : data(rows, cols) {


    threads = calculateThreadsBasedOnMatrixSize(rows,cols, MatrixWrapper::threadScalingFactor,2500);

    Eigen::setNbThreads(threads);

    data.setConstant(value);
}

MatrixWrapper::MatrixWrapper(const Eigen::MatrixXd& other) : data(other) {
     threads = calculateThreadsBasedOnMatrixSize(other.rows(),other.cols(), MatrixWrapper::threadScalingFactor,2500);

    Eigen::setNbThreads(threads);

}


MatrixWrapper::MatrixWrapper(const std::vector<std::vector<double>>& inputData) {
    int rows = inputData.size();
    int cols = rows > 0 ? inputData[0].size() : 0;

    data.resize(rows, cols);

    // Check if all rows have the same number of columns
    for (int i = 0; i < rows; ++i) {
        if (inputData[i].size() != cols) {
            throw std::invalid_argument("All rows must have the same number of columns");
        }
    }

    // Set the number of threads based on hardware concurrency
    // Eigen::setNbThreads(std::thread::hardware_concurrency());

    // Define a lambda function to initialize a range of elements
    auto initRange = [&](int start, int end) {
        for (int i = start; i < end; ++i) {
            data.row(i) = Eigen::Map<const Eigen::RowVectorXd>(&inputData[i][0], cols);
        }
    };

    // Initialize the matrix in parallel using Eigen's internal parallelization
    //#pragma omp parallel for
    for (int i = 0; i < rows; ++i) {
        data.row(i) = Eigen::Map<const Eigen::RowVectorXd>(&inputData[i][0], cols);
    }
}



void addRange(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, Eigen::MatrixXd& C, int startRow, int endRow) {
    C.block(startRow, 0, endRow - startRow, A.cols()) = A.block(startRow, 0, endRow - startRow, A.cols()) + B.block(startRow, 0, endRow - startRow, B.cols());
}

void subtractRange(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, Eigen::MatrixXd& C, int startRow, int endRow) {
    C.block(startRow, 0, endRow - startRow, A.cols()) = A.block(startRow, 0, endRow - startRow, A.cols()) - B.block(startRow, 0, endRow - startRow, B.cols());
}

MatrixWrapper MatrixWrapper::add(const MatrixWrapper& other) const {
    if (data.rows() != other.data.rows() || data.cols() != other.data.cols()) {
        throw std::invalid_argument("Matrix dimensions must match for addition");
    }

    // Enable Eigen's multi-threading
    // Eigen::setNbThreads(std::thread::hardware_concurrency());

    MatrixWrapper result(data.rows(), data.cols());
    result.data = data + other.data;
    return result;
}

MatrixWrapper MatrixWrapper::subtract(const MatrixWrapper& other) const {
    if (data.rows() != other.data.rows() || data.cols() != other.data.cols()) {
        throw std::invalid_argument("Matrix dimensions must match for subtraction");
    }


    // Matrix result(data.rows(), data.cols());
    // result.data = data - other.data;
    // return result;


    return MatrixWrapper(data - other.data);

}


MatrixWrapper MatrixWrapper::log() const {
    MatrixWrapper result(data.rows(), data.cols());
    result.data = data.array().log();
    return result;
}

MatrixWrapper MatrixWrapper::sqrt() const {
    MatrixWrapper result(data.rows(), data.cols());
    result.data = data.array().sqrt();
    return result;
}

MatrixWrapper MatrixWrapper::exp(double base) const {
    MatrixWrapper result(data.rows(), data.cols());
    result.data = data.array().pow(base);
    return result;
}

MatrixWrapper MatrixWrapper::sum(int axis) const {
    if (axis == -1) {
        // Sum all elements and return a 1x1 matrix
        double total_sum = data.sum();
        return MatrixWrapper(1, 1, total_sum);
    } else if (axis == 0) {
        // Sum along columns (resulting in a row vector)
        Eigen::RowVectorXd col_sum = data.colwise().sum();
        MatrixWrapper result(1, col_sum.size());
        result.data = col_sum;
        return result;
    } else if (axis == 1) {
        // Sum along rows (resulting in a column vector)
        Eigen::VectorXd row_sum = data.rowwise().sum();
        MatrixWrapper result(row_sum.size(), 1);
        result.data = row_sum;
        return result;
    } else {
        throw std::invalid_argument("Invalid axis for sum. Axis must be -1, 0, or 1.");
    }
}


MatrixWrapper MatrixWrapper::inverse() const {
    if (data.rows() != data.cols()) {
        throw std::invalid_argument("Matrix must be square to calculate its inverse.");
    }
    return MatrixWrapper(data.inverse().eval());
}

double MatrixWrapper::determinant() const {
    if (data.rows() != data.cols()) {
        throw std::invalid_argument("Matrix must be square to calculate its determinant.");
    }
    return data.determinant();
}

std::pair<MatrixWrapper, MatrixWrapper> MatrixWrapper::eigen() const {
    if (data.rows() != data.cols()) {
        throw std::invalid_argument("Matrix must be square to calculate its eigenvalues and eigenvectors.");
    }
    Eigen::EigenSolver<Eigen::MatrixXd> solver(data);
    MatrixWrapper eigenvalues(1, data.rows());
    MatrixWrapper eigenvectors(data.rows(), data.cols());
    eigenvalues.data = solver.eigenvalues().real();
    eigenvectors.data = solver.eigenvectors().real();
    return std::make_pair(eigenvalues, eigenvectors);
}


MatrixWrapper MatrixWrapper::operator-(const MatrixWrapper& other) const {
    if (data.rows() != other.data.rows() || data.cols() != other.data.cols()) {
        throw std::invalid_argument("Matrix dimensions must match for subtraction");
    }

    return MatrixWrapper(data - other.data);
}


MatrixWrapper MatrixWrapper::operator-(double scalar) const {
    
    return MatrixWrapper(data.array() - scalar);
}

MatrixWrapper MatrixWrapper::transpose() const {
    return MatrixWrapper(data.transpose());
}


void dotProductRange(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, Eigen::MatrixXd& C, int startRow, int endRow) {
    C.block(startRow, 0, endRow - startRow, B.cols()) = A.block(startRow, 0, endRow - startRow, A.cols()) * B;
}


void dotProductRange(const Eigen::Ref<const Eigen::MatrixXd>& data,
                                const Eigen::Ref<const Eigen::MatrixXd>& otherData,
                                Eigen::Ref<Eigen::MatrixXd> resultData,
                                int startRow, int endRow) {
        resultData.block(startRow, 0, endRow - startRow, resultData.cols()) =
            data.block(startRow, 0, endRow - startRow, data.cols()) * otherData;
    }


MatrixWrapper MatrixWrapper::dot(const MatrixWrapper& other) const {
    if (data.cols() != other.data.rows()) {
        throw std::invalid_argument("Matrix dimensions must be compatible for dot product");
    }

    MatrixWrapper result(data.rows(), other.data.cols());
    int numThreads = threads;
    int rowsPerThread = data.rows() / numThreads;

    ThreadManager threadManager(threads);

    auto dataRef = std::cref(data);
    auto otherDataRef = std::cref(other.data);
    auto resultDataRef = std::ref(result.data);

    for (int i = 0; i < numThreads; ++i) {
        int startRow = i * rowsPerThread;
        int endRow = (i == numThreads - 1) ? data.rows() : (i + 1) * rowsPerThread;

        threadManager.submitTask([startRow, endRow, dataRef, otherDataRef, resultDataRef] {
            dotProductRange(dataRef, otherDataRef, resultDataRef, startRow, endRow);
        });
    }

    threadManager.waitForCompletion();

    return result;
}

// Matrix Matrix::dot(const Matrix& other) const {
//     if (data.cols() != other.data.rows()) {
//         throw std::invalid_argument("Matrix dimensions must be compatible for dot product");
//     }

//     Matrix result(data.rows(), other.data.cols());
//     int numThreads = threads;
//     int rowsPerThread = data.rows() / numThreads;

//     auto dataRef = std::cref(data);
//     auto otherDataRef = std::cref(other.data);
//     auto resultDataRef = std::ref(result.data);

//     std::vector<std::thread> threads;
//     for (int i = 0; i < numThreads; ++i) {
//         int startRow = i * rowsPerThread;
//         int endRow = (i == numThreads - 1) ? data.rows() : (i + 1) * rowsPerThread;
//         threads.push_back(std::thread(dotProductRange, dataRef, otherDataRef, resultDataRef, startRow, endRow));
//     }

//     for (auto& t : threads) {
//         t.join();
//     }

//     return result;
// }

MatrixWrapper MatrixWrapper::addScalar(const double other) const {
    MatrixWrapper result(data.rows(), data.cols());
    result.data = data.array() + other;
    return result;
}

MatrixWrapper MatrixWrapper::subtractScalar(const double other) const {
    MatrixWrapper result(data.rows(), data.cols());
    result.data = data.array() - other;
    return result;
}


// other functions 
std::vector<int> MatrixWrapper::argmax(int axis) const {
    std::vector<int> indices;

    if (axis == 0) {
        // Max along columns
        indices.resize(data.cols());
        // Eigen::setNbThreads(this->threads);
        
        #pragma omp parallel for num_threads(this->threads)
        for (int j = 0; j < data.cols(); ++j) {
            int maxIndex = 0;
            for (int i = 1; i < data.rows(); ++i) {
                if (data(i, j) > data(maxIndex, j)) {
                    maxIndex = i;
                }
            }
            indices[j] = maxIndex;
        }

    } else if (axis == 1) {
        // Max along rows
        indices.resize(data.rows());
        Eigen::setNbThreads(this->threads);
        
        #pragma omp parallel for num_threads(this->threads)
        for (int i = 0; i < data.rows(); ++i) {
            int maxIndex = 0;
            for (int j = 1; j < data.cols(); ++j) {
                if (data(i, j) > data(i, maxIndex)) {
                    maxIndex = j;
                }
            }
            indices[i] = maxIndex;
        }

    } else {
        throw std::invalid_argument("Axis must be 0 or 1");
    }

    return indices;
}


MatrixWrapper MatrixWrapper::clip(double min_val, double max_val) const {
    MatrixWrapper result(data.rows(), data.cols());
    result.data = data.array().min(max_val).max(min_val);
    return result;
}


    // Operator overloads for +, -, *

MatrixWrapper MatrixWrapper::operator/(const MatrixWrapper& other) const {
    if (data.rows() != other.getRows() || data.cols() != other.getCols()) {
        throw std::invalid_argument("Matrices must have the same dimensions for division.");
    }
    
    return MatrixWrapper(data.array() / other.data.array());
}

MatrixWrapper MatrixWrapper::operator/(double scalar) const {
    if (scalar == 0) {
        throw std::invalid_argument("Division by zero encountered in scalar division.");
    }

    return MatrixWrapper(data.array() + scalar);

}


MatrixWrapper MatrixWrapper::operator*(const MatrixWrapper& other) const {
    if (data.rows() != other.getRows() || data.cols() != other.getCols()) {
        throw std::invalid_argument("Matrices must have the same dimensions for division.");
    }
    
    return MatrixWrapper(data.array() * other.data.array());
}

MatrixWrapper MatrixWrapper::operator*(double scalar) const {
    if (scalar == 0) {
        throw std::invalid_argument("Division by zero encountered in scalar division.");
    }

    return MatrixWrapper(data.array() * scalar);

}


MatrixWrapper MatrixWrapper::random(int rows, int cols, double min, double max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);

    Eigen::MatrixXd matrix = Eigen::MatrixXd::NullaryExpr(rows, cols, [&]() { return dis(gen); });

    return MatrixWrapper(matrix);
}


void MatrixWrapper::display() const {
    std::cout << data << std::endl;
}



// Helper functions 

// Function to calculate a safe number of threads for context switching based on matrix size
int calculateThreadsBasedOnMatrixSize(int rows, int cols, double scalingFactor = 1.0, int maxThreads = 2500) {
    // Calculate the total number of elements in the matrix
    int totalElements = rows * cols;

    // Apply logarithm base 2 transformation
    double logElements = std::log2(static_cast<double>(totalElements));

    // Calculate the number of threads
    int numThreads = static_cast<int>(logElements * scalingFactor);

    // Get the number of hardware threads
    int hardwareThreads = std::thread::hardware_concurrency();
    if (hardwareThreads == 0) {
        // If hardware_concurrency() returns 0, assume a default value of 4
        hardwareThreads = 4;
    }

    // Ensure the number of threads is at least the number of hardware threads
    numThreads = std::max(numThreads, hardwareThreads);

    // Cap the number of threads to a maximum limit
    return std::min(numThreads, maxThreads);
}
