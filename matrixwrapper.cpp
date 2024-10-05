#include "matrixwrapper.h"
#include <stdexcept>
#include <iostream>
#include <thread>
#include <vector>
#include <cmath> 
#include <functional>
#include "ThreadManager.h"
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

double MatrixWrapper::getItem(int row, int col) const {
    if (row < 0 || row >= data.rows() || col < 0 || col >= data.cols()) {
        throw std::out_of_range("Matrix indices out of range");
    }
    return data(row, col);
}

void MatrixWrapper::setItem(int row, int col, double value) {
    if (row < 0 || row >= data.rows() || col < 0 || col >= data.cols()) {
        throw std::out_of_range("Matrix indices out of range");
    }
    data(row, col) = value;
}

MatrixWrapper MatrixWrapper::getRow(int row) const {
    if (row < 0 || row >= data.rows()) {
        throw std::out_of_range("Row index out of range");
    }
    return MatrixWrapper(data.row(row));
}

void MatrixWrapper::setRow(int row, const MatrixWrapper& rowData) {
    if (row < 0 || row >= data.rows()) {
        throw std::out_of_range("Row index out of range");
    }
    if (rowData.data.cols() != data.cols()) {
        throw std::invalid_argument("Row data must have the same number of columns as the matrix");
    }
    data.row(row) = rowData.data;
}

MatrixWrapper MatrixWrapper::getCol(int col) const {
    if (col < 0 || col >= data.cols()) {
        throw std::out_of_range("Column index out of range");
    }
    return MatrixWrapper(data.col(col));
}

void MatrixWrapper::setCol(int col, const MatrixWrapper& colData) {
    if (col < 0 || col >= data.cols()) {
        throw std::out_of_range("Column index out of range");
    }
    if (colData.data.rows() != data.rows()) {
        throw std::invalid_argument("Column data must have the same number of rows as the matrix");
    }
    data.col(col) = colData.data;
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
    // Case 1: Dimensions match exactly
    if (data.rows() == other.data.rows() && data.cols() == other.data.cols()) {
        MatrixWrapper result(data.rows(), data.cols());
        result.data = data + other.data;
        return result;
    }
    // Case 2: This matrix has 1 row, other has multiple rows
    else if (data.rows() == 1 && data.cols() == other.data.cols()) {
        MatrixWrapper result(other.data.rows(), data.cols());
        result.data = other.data.rowwise() + data.row(0);
        return result;
    }
    // Case 3: Other matrix has 1 row, this has multiple rows
    else if (other.data.rows() == 1 && data.cols() == other.data.cols()) {
        MatrixWrapper result(data.rows(), data.cols());
        result.data = data.rowwise() + other.data.row(0);
        return result;
    }
    // Case 4: This matrix has 1 column, other has multiple columns
    else if (data.cols() == 1 && data.rows() == other.data.rows()) {
        MatrixWrapper result(data.rows(), other.data.cols());
        result.data = other.data.colwise() + data.col(0);
        return result;
    }
    // Case 5: Other matrix has 1 column, this has multiple columns
    else if (other.data.cols() == 1 && other.data.rows() == data.rows()) {
        MatrixWrapper result(data.rows(), data.cols());
        result.data = data.colwise() + other.data.col(0);
        return result;
    }
    else {
        throw std::invalid_argument("Matrix dimensions are not compatible for addition with broadcasting");
    }
}

MatrixWrapper MatrixWrapper::subtract(const MatrixWrapper &other) const {
    if (data.rows() == other.data.rows() && data.cols() == other.data.cols()) {
        return MatrixWrapper(data - other.data);
    } else if (data.cols() == other.data.cols() && other.data.rows() == 1) {
        return MatrixWrapper(data.rowwise() - other.data.row(0));
    } else if (data.rows() == other.data.rows() && other.data.cols() == 1) {
        return MatrixWrapper(data.colwise() - other.data.col(0));
    } else {
        throw std::invalid_argument("Matrix dimensions do not match for broadcasting");
    }
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
    result.data = data.array().exp();
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
    for (int i = 0; i < data.rows(); ++i) {
        for (int j = 0; j < data.cols(); ++j) {
            result.data(i, j) = std::min(std::max(data(i, j), min_val), max_val);
        }
    }
    return result;
}


    // Operator overloads for +, -, *

// MatrixWrapper MatrixWrapper::operator/(const MatrixWrapper& other) const {
//     if (data.rows() != other.getRows() || data.cols() != other.getCols()) {
//         throw std::invalid_argument("Matrices must have the same dimensions for division.");
//     }
    
//     return MatrixWrapper(data.array() / other.data.array());
// }

MatrixWrapper MatrixWrapper::operator/(const MatrixWrapper& other) const {
    if (data.cols() != other.data.cols() && other.data.cols() != 1) {
        throw std::invalid_argument("Number of columns must match for division or divisor must be a column vector");
    }

    MatrixWrapper result(data.rows(), data.cols());

    if (data.rows() == other.data.rows() && data.cols() == other.data.cols()) {
        // Element-wise division
        result.data = data.array() / other.data.array();
    } else if (other.data.rows() == 1 && other.data.cols() == data.cols()) {
        // Broadcasting across rows
        for (int i = 0; i < data.rows(); ++i) {
            result.data.row(i) = data.row(i).array() / other.data.row(0).array();
        }
    } else if (other.data.cols() == 1 && other.data.rows() == data.rows()) {
        // Broadcasting across columns
        for (int i = 0; i < data.cols(); ++i) {
            result.data.col(i) = data.col(i).array() / other.data.col(0).array();
        }
    } else {
        throw std::invalid_argument("Matrices must have compatible dimensions for division with broadcasting");
    }

    return result;
}

MatrixWrapper MatrixWrapper::operator/(double scalar) const {
    if (scalar == 0) {
        throw std::invalid_argument("Division by zero encountered in scalar division.");
    }

    return MatrixWrapper(data.array() / scalar);

}


MatrixWrapper MatrixWrapper::operator*(const MatrixWrapper& other) const {
    if (data.rows() != other.getRows() || data.cols() != other.getCols()) {
        throw std::invalid_argument("Matrices must have the same dimensions for division.");
    }
    
    return MatrixWrapper(data.array() * other.data.array());
}

MatrixWrapper MatrixWrapper::operator*(double scalar) const {

    return MatrixWrapper(data.array() * scalar);

}

MatrixWrapper MatrixWrapper::random(int rows, int cols, double min, double max) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);

    Eigen::MatrixXd matrix = Eigen::MatrixXd::NullaryExpr(rows, cols, [&]() { return dis(gen); });

    return MatrixWrapper(matrix);
}

MatrixWrapper MatrixWrapper::randomBinomial(int rows, int cols) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::bernoulli_distribution dist;

    Eigen::MatrixXd matrix(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix(i, j) = dist(gen) ? 1.0 : 0.0;
        }
    }

    return MatrixWrapper(matrix);
}

MatrixWrapper MatrixWrapper::relu(double condition, double val) const {
    MatrixWrapper result(data.rows(), data.cols());
    result.data = data.unaryExpr([condition, val](double x) {
        return x < condition ? val : x;
    });
    return result;
}
// MatrixWrapper MatrixWrapper::relu() const {
//     MatrixWrapper result(data.rows(), data.cols());
//     Eigen::setNbThreads(threads); // Set the number of threads
//     result.data = data.parallelUnaryExpr([](double x) {
//         return std::max(x, 0.0);
//     });
//     return result;
// }

MatrixWrapper MatrixWrapper::multiply(const MatrixWrapper &other) const {
     if (data.rows() != other.data.rows() || data.cols() != other.data.cols()) {
        throw std::invalid_argument("Matrix dimensions must match for element-wise multiplication.");
    }
    return MatrixWrapper(data.cwiseProduct(other.data));
}

MatrixWrapper MatrixWrapper::multiplyElementwise(const MatrixWrapper &other) const {
    if (data.rows() != other.data.rows() && data.cols() != other.data.cols() && other.data.cols() != 1) {
        throw std::invalid_argument("Matrix dimensions must agree for element-wise multiplication or allow for broadcasting.");
    }
    
    MatrixWrapper result(data.rows(), data.cols());
    
    if (other.data.cols() == 1) {
        // Broadcasting case
        for (int i = 0; i < data.rows(); ++i) {
            for (int j = 0; j < data.cols(); ++j) {
                result.data(i, j) = data(i, j) * other.data(i, 0);
            }
        }
    } else {
        // Standard element-wise multiplication
        result.data = data.array() * other.data.array();
    }
    
    return result;
}

MatrixWrapper MatrixWrapper::std(int axis, double ddof) const {
    Eigen::MatrixXd centered;
    Eigen::MatrixXd variance;

    if (axis == 0) {
        Eigen::RowVectorXd mean = data.colwise().mean();
        centered = data.rowwise() - mean;
        variance = centered.array().square().colwise().sum() / (data.rows() - ddof);
        return MatrixWrapper(variance.array().sqrt().matrix());
    } else if (axis == 1) {
        Eigen::VectorXd mean = data.rowwise().mean();
        centered = data.colwise() - mean;
        variance = centered.array().square().rowwise().sum() / (data.cols() - ddof);
        return MatrixWrapper(variance.array().sqrt().matrix());
    } else {
        double mean = data.mean();
        centered = data.array() - mean;
        double var = centered.array().square().sum() / (data.size() - ddof);
        return MatrixWrapper(Eigen::MatrixXd::Constant(1, 1, std::sqrt(var)));
    }
}



void MatrixWrapper::display() const {
    std::cout << data << std::endl;
}

MatrixWrapper MatrixWrapper::reshape(const std::vector<int>& newshape, const char order) const {
    int new_rows = newshape.size() > 0 ? newshape[0] : 1;
    int new_cols = data.size() / new_rows;

    if (new_rows * new_cols != data.size()) {
        throw std::invalid_argument("The new shape is not compatible with the current shape and number of elements.");
    }

    MatrixWrapper result(new_rows, new_cols);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> result_map(result.data.data(), new_rows, new_cols);

    if (order == 'C') {
        for (int i = 0; i < data.size(); ++i) {
            result_map(i) = data(i);
        }
    } else if (order == 'F') {
        for (int i = 0; i < new_rows; ++i) {
            for (int j = 0; j < new_cols; ++j) {
                result_map(i, j) = data(j * new_rows + i);
            }
        }
    } else {
        throw std::invalid_argument("Invalid order. Order must be 'C' or 'F'.");
    }

    return result;
}


MatrixWrapper MatrixWrapper::abs() const {
    MatrixWrapper result(data.rows(), data.cols());
    result.data = data.array().abs();
    return result;
}


MatrixWrapper MatrixWrapper::sum(int axis, double initial, const Eigen::MatrixXd& where) const {
    if (axis == -1) {
        // Sum all elements and return a 1x1 matrix
        double total_sum = data.sum();
        if (!where.isZero()) {
            total_sum = (data.array() * where.array()).sum();
        }
        total_sum += initial;
        return MatrixWrapper(1, 1, total_sum);
    } else if (axis == 0) {
        // Sum along columns (resulting in a row vector)
        Eigen::RowVectorXd col_sum = data.colwise().sum();
        if (!where.isZero()) {
            col_sum = (data.array() * where.array()).colwise().sum();
        }
        col_sum.array() += initial;
        MatrixWrapper result(1, col_sum.size());
        result.data = col_sum;
        return result;
    } else if (axis == 1) {
        // Sum along rows (resulting in a column vector)
        Eigen::VectorXd row_sum = data.rowwise().sum();
        if (!where.isZero()) {
            row_sum = (data.array() * where.array()).rowwise().sum();
        }
        row_sum.array() += initial;
        MatrixWrapper result(row_sum.size(), 1);
        result.data = row_sum;
        return result;
    } else {
        throw std::invalid_argument("Invalid axis for sum. Axis must be -1, 0, or 1.");
    }
}


MatrixWrapper MatrixWrapper::pow(const MatrixWrapper& exponent) const {
    if (data.rows() != exponent.data.rows() || data.cols() != exponent.data.cols()) {
        throw std::invalid_argument("Base and exponent matrices must have the same dimensions");
    }

    MatrixWrapper result(data.rows(), data.cols());
    result.data = data.array().pow(exponent.data.array());
    return result;
}

MatrixWrapper MatrixWrapper::pow(double exponent) const {
    MatrixWrapper result(data.rows(), data.cols());
    result.data = data.array().pow(exponent);
    return result;
}

MatrixWrapper MatrixWrapper::getValuesFromIndices(const std::vector<int>& indices) const {

    std::vector<std::vector<double>> resultData;

    std::vector<double> values;

    int row = 0;

    for (int index : indices) {

        if (row < data.rows() && index >= 0 && index < data.cols()) {
            values.push_back(data(row, index));
        } else {
            values.push_back(std::nan("")); // Use std::nan("") for null value
        }

        row++;
    }

    resultData.push_back(values);

    return MatrixWrapper(resultData);
}


MatrixWrapper MatrixWrapper::eye(int N, int M, int k, const std::string& dtype, const std::string& order) {
    if (M == -1) {
        M = N;
    }

    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(N, M);

    if (k >= 0) {
        for (int i = 0; i < std::min(N, M - k); ++i) {
            matrix(i, i + k) = 1.0;
        }
    } else {
        for (int i = 0; i < std::min(N + k, M); ++i) {
            matrix(i - k, i) = 1.0;
        }
    }

    return MatrixWrapper(matrix);
}


MatrixWrapper MatrixWrapper::selectRowsByIndices(const std::vector<int>& indices) const {
    int numRows = indices.size();
    int numCols = data.cols();
    Eigen::MatrixXd resultData(numRows, numCols);

    for (int i = 0; i < numRows; ++i) {
        int index = indices[i];
        if (index >= 0 && index < data.rows()) {
            resultData.row(i) = data.row(index);
        } else {
            resultData.row(i).setConstant(std::nan(""));
        }
    }

    return MatrixWrapper(resultData);
}

MatrixWrapper MatrixWrapper::mean(int axis, const std::string& dtype, const Eigen::MatrixXd& where) const {
    if (axis == -1) {
        // Calculate the mean of all elements
        Eigen::MatrixXd matrixToMean = where.isZero() ? data : data.cwiseProduct(where);
        double totalSum = matrixToMean.sum();
        int count = where.isZero() ? matrixToMean.size() : where.count();
        return MatrixWrapper(1, 1, totalSum / count);
    } else if (axis == 0) {
        // Calculate the mean along columns
        Eigen::MatrixXd matrixToMean = where.isZero() ? data : data.cwiseProduct(where);
        Eigen::RowVectorXd colMeans = Eigen::RowVectorXd::Zero(matrixToMean.cols());
        for (int j = 0; j < matrixToMean.cols(); ++j) {
            colMeans(j) = matrixToMean.col(j).sum() / matrixToMean.col(j).count();
        }
        MatrixWrapper result(1, colMeans.size());
        result.data = colMeans;
        return result;
    } else if (axis == 1) {
        // Calculate the mean along rows
        Eigen::MatrixXd matrixToMean = where.isZero() ? data : data.cwiseProduct(where);
        Eigen::VectorXd rowMeans = Eigen::VectorXd::Zero(matrixToMean.rows());
        for (int i = 0; i < matrixToMean.rows(); ++i) {
            rowMeans(i) = matrixToMean.row(i).sum() / matrixToMean.row(i).count();
        }
        MatrixWrapper result(rowMeans.size(), 1);
        result.data = rowMeans;
        return result;
    } else {
        throw std::invalid_argument("Invalid axis for mean. Axis must be -1, 0, or 1.");
    }
}

MatrixWrapper MatrixWrapper::sign(const Eigen::MatrixXd& where) const {
    MatrixWrapper result(data.rows(), data.cols());
    if (where.size() == 0 || where.isOnes()) {
        result.data = data.unaryExpr([](double x) { return (x > 0 ? 1.0 : (x < 0 ? -1.0 : 0.0)); });
    } else {
        result.data = (data.array() > 0).select(where.cast<double>(), 0.0) - (data.array() < 0).select(where.cast<double>(), 0.0);
    }
    return result;
}


MatrixWrapper MatrixWrapper::copy() const {
    MatrixWrapper newMatrix(this->data.rows(), this->data.cols());
    newMatrix.data = this->data;  // This performs a deep copy of the Eigen matrix
    return newMatrix;
}


MatrixWrapper MatrixWrapper::max(int axis, double initial, const Eigen::MatrixXd& where) const {
    if (axis == -1) {
        // Maximum of all elements
        double maxVal = where.isZero() ? data.maxCoeff() : (data.array() * where.array()).maxCoeff();
        maxVal = std::max(maxVal, initial);
        return MatrixWrapper(1, 1, maxVal);
    } else if (axis == 0) {
        // Maximum along columns
        Eigen::RowVectorXd maxVals = data.colwise().maxCoeff();
        if (!where.isZero()) {
            maxVals = (data.array() * where.array()).colwise().maxCoeff();
        }
        for (int i = 0; i < maxVals.size(); ++i) {
            maxVals(i) = std::max(maxVals(i), initial);
        }
        MatrixWrapper result(1, maxVals.size());
        result.data = maxVals;
        return result;
    } else if (axis == 1) {
        // Maximum along rows
        Eigen::VectorXd maxVals = data.rowwise().maxCoeff();
        if (!where.isZero()) {
            maxVals = (data.array() * where.array()).rowwise().maxCoeff();
        }
        for (int i = 0; i < maxVals.size(); ++i) {
            maxVals(i) = std::max(maxVals(i), initial);
        }
        MatrixWrapper result(maxVals.size(), 1);
        result.data = maxVals;
        return result;
    } else {
        throw std::invalid_argument("Invalid axis for max. Axis must be -1, 0, or 1.");
    }
}


MatrixWrapper MatrixWrapper::min(int axis, double initial, const Eigen::MatrixXd& where) const {
    if (axis == -1) {
        // Minimum of all elements
        double minVal = where.isZero() ? data.minCoeff() : (data.array() * where.array()).minCoeff();
        minVal = std::min(minVal, initial);
        return MatrixWrapper(1, 1, minVal);
    } else if (axis == 0) {
        // Minimum along columns
        Eigen::RowVectorXd minVals = data.colwise().minCoeff();
        if (!where.isZero()) {
            minVals = (data.array() * where.array()).colwise().minCoeff();
        }
        for (int i = 0; i < minVals.size(); ++i) {
            minVals(i) = std::min(minVals(i), initial);
        }
        MatrixWrapper result(1, minVals.size());
        result.data = minVals;
        return result;
    } else if (axis == 1) {
        // Minimum along rows
        Eigen::VectorXd minVals = data.rowwise().minCoeff();
        if (!where.isZero()) {
            minVals = (data.array() * where.array()).rowwise().minCoeff();
        }
        for (int i = 0; i < minVals.size(); ++i) {
            minVals(i) = std::min(minVals(i), initial);
        }
        MatrixWrapper result(minVals.size(), 1);
        result.data = minVals;
        return result;
    } else {
        throw std::invalid_argument("Invalid axis for min. Axis must be -1, 0, or 1.");
    }
}

MatrixWrapper MatrixWrapper::glorot_uniform(int fan_in, int fan_out) {
    double limit = std::sqrt(6.0 / (fan_in + fan_out));

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-limit, limit);

    MatrixWrapper result(fan_in, fan_out);
    for (int i = 0; i < fan_in; ++i) {
        for (int j = 0; j < fan_out; ++j) {
            result.data(i, j) = dis(gen);
        }
    }

    return result;
}

// In MatrixWrapper.cpp

MatrixWrapper MatrixWrapper::getSlice(const std::vector<int>& range) const {
    if (range.size() != 2) {
        throw std::invalid_argument("Slice range must contain exactly two elements");
    }
    int start = range[0];
    int end = range[1];
    if (start < 0 || end > data.rows() || start >= end) {
        throw std::out_of_range("Invalid slice range");
    }
    return MatrixWrapper(data.block(start, 0, end - start, data.cols()));
}

MatrixWrapper MatrixWrapper::zeros_like(const MatrixWrapper& other) {
    return MatrixWrapper(other.getRows(), other.getCols(), 0.0);
}

MatrixWrapper MatrixWrapper::zeros_like() const {
    return MatrixWrapper(this->getRows(), this->getCols(), 0.0);
}

MatrixWrapper MatrixWrapper::ones_like(const MatrixWrapper& other) {
    return MatrixWrapper(other.getRows(), other.getCols(), 1.0);
}

MatrixWrapper MatrixWrapper::ones_like() const {
    return MatrixWrapper(this->getRows(), this->getCols(), 1.0);
}

MatrixWrapper MatrixWrapper::round(int precision) {
    double scale = std::pow(10.0, precision);
    MatrixWrapper result(data.rows(), data.cols());
    result.data = (data.array() * scale).round() / scale;
    return result;
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
