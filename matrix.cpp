#include "matrix.h"
#include <stdexcept>
#include <iostream>
#include <thread>
#include <vector>
#include <cmath> 
#include <functional>

Matrix::Matrix(int rows, int cols, double value) : data(rows, cols) {
    data.setConstant(value);
}

Matrix::Matrix(const Eigen::MatrixXd& other) : data(other) {}

Matrix::Matrix(const std::vector<std::vector<double>>& inputData) {
    int rows = inputData.size();
    int cols = rows > 0 ? inputData[0].size() : 0;
    data.resize(rows, cols);

    // Check if all rows have the same number of columns
    for (int i = 0; i < rows; ++i) {
        if (inputData[i].size() != cols) {
            throw std::invalid_argument("All rows must have the same number of columns");
        }
    }

    // Log the number of rows and columns
    std::cout << "Initializing matrix with " << rows << " rows and " << cols << " columns using multi-threading." << std::endl;

    // Determine the number of threads
    int numThreads = std::min(20, static_cast<int>(std::thread::hardware_concurrency()));
    int rowsPerThread = rows / numThreads;

    // Define a lambda function to initialize a range of elements
    auto initRange = [&](int startRow, int endRow) {
        for (int i = startRow; i < endRow; ++i) {
            for (int j = 0; j < cols; ++j) {
                data(i, j) = inputData[i][j];
            }
        }
        // Debug: Log each thread's range
        std::cout << "Thread completed rows " << startRow << " to " << endRow << std::endl;
    };

    // Create and run threads
    std::vector<std::thread> threads;
    for (int i = 0; i < numThreads; ++i) {
        int startRow = i * rowsPerThread;
        int endRow = (i == numThreads - 1) ? rows : (i + 1) * rowsPerThread;
        threads.emplace_back(initRange, startRow, endRow);
    }

    // Wait for all threads to finish
    for (auto& thread : threads) {
        if (thread.joinable()) {
            thread.join();
        }
    }

    std::cout << "Matrix initialization complete." << std::endl;
}


void addRange(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, Eigen::MatrixXd& C, int startRow, int endRow) {
    C.block(startRow, 0, endRow - startRow, A.cols()) = A.block(startRow, 0, endRow - startRow, A.cols()) + B.block(startRow, 0, endRow - startRow, B.cols());
}

void subtractRange(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, Eigen::MatrixXd& C, int startRow, int endRow) {
    C.block(startRow, 0, endRow - startRow, A.cols()) = A.block(startRow, 0, endRow - startRow, A.cols()) - B.block(startRow, 0, endRow - startRow, B.cols());
}

Matrix Matrix::add(const Matrix& other) const {
    if (data.rows() != other.data.rows() || data.cols() != other.data.cols()) {
        throw std::invalid_argument("Matrix dimensions must match for addition");
    }

    // Enable Eigen's multi-threading
    Eigen::setNbThreads(std::thread::hardware_concurrency());

    Matrix result(data.rows(), data.cols());
    result.data = data + other.data;
    return result;
}

Matrix Matrix::subtract(const Matrix& other) const {
    if (data.rows() != other.data.rows() || data.cols() != other.data.cols()) {
        throw std::invalid_argument("Matrix dimensions must match for subtraction");
    }

    // Enable Eigen's multi-threading
    Eigen::setNbThreads(std::thread::hardware_concurrency());

    Matrix result(data.rows(), data.cols());
    result.data = data - other.data;
    return result;
}


Matrix Matrix::log() const {
    Matrix result(data.rows(), data.cols());
    result.data = data.array().log();
    return result;
}

Matrix Matrix::sqrt() const {
    Matrix result(data.rows(), data.cols());
    result.data = data.array().sqrt();
    return result;
}

Matrix Matrix::exp(double base) const {
    Matrix result(data.rows(), data.cols());
    result.data = data.array().pow(base);
    return result;
}

Matrix Matrix::sum(int axis) const {
    if (axis == -1) {
        // Sum all elements and return a 1x1 matrix
        double total_sum = data.sum();
        return Matrix(1, 1, total_sum);
    } else if (axis == 0) {
        // Sum along columns (resulting in a row vector)
        Eigen::RowVectorXd col_sum = data.colwise().sum();
        Matrix result(1, col_sum.size());
        result.data = col_sum;
        return result;
    } else if (axis == 1) {
        // Sum along rows (resulting in a column vector)
        Eigen::VectorXd row_sum = data.rowwise().sum();
        Matrix result(row_sum.size(), 1);
        result.data = row_sum;
        return result;
    } else {
        throw std::invalid_argument("Invalid axis for sum. Axis must be -1, 0, or 1.");
    }
}


Matrix Matrix::inverse() const {
    if (data.rows() != data.cols()) {
        throw std::invalid_argument("Matrix must be square to calculate its inverse.");
    }
    return Matrix(data.inverse().eval());
}

double Matrix::determinant() const {
    if (data.rows() != data.cols()) {
        throw std::invalid_argument("Matrix must be square to calculate its determinant.");
    }
    return data.determinant();
}

std::pair<Matrix, Matrix> Matrix::eigen() const {
    if (data.rows() != data.cols()) {
        throw std::invalid_argument("Matrix must be square to calculate its eigenvalues and eigenvectors.");
    }
    Eigen::EigenSolver<Eigen::MatrixXd> solver(data);
    Matrix eigenvalues(1, data.rows());
    Matrix eigenvectors(data.rows(), data.cols());
    eigenvalues.data = solver.eigenvalues().real();
    eigenvectors.data = solver.eigenvectors().real();
    return std::make_pair(eigenvalues, eigenvectors);
}


Matrix Matrix::operator*(double scalar) const {
    Matrix result(data.rows(), data.cols());
    result.data = data * scalar;
    return result;
}


Matrix Matrix::transpose() const {
    Matrix result(data.cols(), data.rows());
    result.data = data.transpose();
    return result;
}


void dotProductRange(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, Eigen::MatrixXd& C, int startRow, int endRow) {
    C.block(startRow, 0, endRow - startRow, B.cols()) = A.block(startRow, 0, endRow - startRow, A.cols()) * B;
}

Matrix Matrix::dot(const Matrix& other) const {
    if (data.cols() != other.data.rows()) {
        throw std::invalid_argument("Matrix dimensions must be compatible for dot product");
    }

    Matrix result(data.rows(), other.data.cols());
    int numThreads = std::min(20, static_cast<int>(std::thread::hardware_concurrency())); // Use at most 20 threads or the number of available hardware threads
    int rowsPerThread = data.rows() / numThreads;

    auto dataRef = std::cref(data);
    auto otherDataRef = std::cref(other.data);
    auto resultDataRef = std::ref(result.data);

    std::vector<std::thread> threads;
    for (int i = 0; i < numThreads; ++i) {
        int startRow = i * rowsPerThread;
        int endRow = (i == numThreads - 1) ? data.rows() : (i + 1) * rowsPerThread;
        threads.push_back(std::thread(dotProductRange, dataRef, otherDataRef, resultDataRef, startRow, endRow));
    }

    for (auto& t : threads) {
        t.join();
    }

    return result;
}

Matrix Matrix::addScalar(const double other) const {
    Matrix result(data.rows(), data.cols());
    result.data = data.array() + other;
    return result;
}

Matrix Matrix::subtractScalar(const double other) const {
    Matrix result(data.rows(), data.cols());
    result.data = data.array() - other;
    return result;
}


// other functions 
std::vector<int> Matrix::argmax(int axis) const {
    std::vector<int> indices;

    if (axis == 0) {
        // Max along columns
        indices.resize(data.cols());

        auto find_max_in_col = [&](int start, int end) {
            for (int j = start; j < end; ++j) {
                int maxIndex = 0;
                for (int i = 1; i < data.rows(); ++i) {
                    if (data(i, j) > data(maxIndex, j)) {
                        maxIndex = i;
                    }
                }
                indices[j] = maxIndex;
            }
        };

        int num_threads = std::thread::hardware_concurrency();
        std::vector<std::thread> threads(num_threads);

        int cols_per_thread = data.cols() / num_threads;
        int remaining_cols = data.cols() % num_threads;

        int start = 0;
        for (int i = 0; i < num_threads; ++i) {
            int end = start + cols_per_thread + (remaining_cols > 0 ? 1 : 0);
            if (remaining_cols > 0) remaining_cols--;
            threads[i] = std::thread(find_max_in_col, start, end);
            start = end;
        }

        for (auto& th : threads) {
            th.join();
        }

    } else if (axis == 1) {
        // Max along rows
        indices.resize(data.rows());

        auto find_max_in_row = [&](int start, int end) {
            for (int i = start; i < end; ++i) {
                int maxIndex = 0;
                for (int j = 1; j < data.cols(); ++j) {
                    if (data(i, j) > data(i, maxIndex)) {
                        maxIndex = j;
                    }
                }
                indices[i] = maxIndex;
            }
        };

        int num_threads = std::thread::hardware_concurrency();
        std::vector<std::thread> threads(num_threads);

        int rows_per_thread = data.rows() / num_threads;
        int remaining_rows = data.rows() % num_threads;

        int start = 0;
        for (int i = 0; i < num_threads; ++i) {
            int end = start + rows_per_thread + (remaining_rows > 0 ? 1 : 0);
            if (remaining_rows > 0) remaining_rows--;
            threads[i] = std::thread(find_max_in_row, start, end);
            start = end;
        }

        for (auto& th : threads) {
            th.join();
        }

    } else {
        throw std::invalid_argument("Axis must be 0 or 1");
    }

    return indices;
}


Matrix Matrix::clip(double min_val, double max_val) const {
    Matrix result(data.rows(), data.cols());
    result.data = data.array().min(max_val).max(min_val);
    return result;
}


    // Operator overloads for +, -, *

Matrix Matrix::operator/(const Matrix& other) const {
    if (data.rows() != other.getRows() || data.cols() != other.getCols()) {
        throw std::invalid_argument("Matrices must have the same dimensions for division.");
    }
    Matrix result(data.rows(), data.cols());

    auto divide_row = [&](int start, int end) {
        for (int i = start; i < end; ++i) {
            for (int j = 0; j < data.cols(); ++j) {
                if (other(i, j) == 0) {
                    throw std::invalid_argument("Division by zero encountered in matrix division.");
                }
                result(i, j) = data(i, j) / other(i, j);
            }
        }
    };

    int num_threads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads(num_threads);

    int rows_per_thread = data.rows() / num_threads;
    int remaining_rows = data.rows() % num_threads;

    int start = 0;
    for (int i = 0; i < num_threads; ++i) {
        int end = start + rows_per_thread + (remaining_rows > 0 ? 1 : 0);
        if (remaining_rows > 0) remaining_rows--;
        threads[i] = std::thread(divide_row, start, end);
        start = end;
    }

    for (auto& th : threads) {
        th.join();
    }

    return result;
}

Matrix Matrix::operator/(double scalar) const {
    if (scalar == 0) {
        throw std::invalid_argument("Division by zero encountered in scalar division.");
    }
    Matrix result(data.rows(), data.cols());

    auto divide_row = [&](int start, int end) {
        for (int i = start; i < end; ++i) {
            for (int j = 0; j < data.cols(); ++j) {
                result(i, j) = data(i, j) / scalar;
            }
        }
    };

    int num_threads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads(num_threads);

    int rows_per_thread = data.rows() / num_threads;
    int remaining_rows = data.rows() % num_threads;

    int start = 0;
    for (int i = 0; i < num_threads; ++i) {
        int end = start + rows_per_thread + (remaining_rows > 0 ? 1 : 0);
        if (remaining_rows > 0) remaining_rows--;
        threads[i] = std::thread(divide_row, start, end);
        start = end;
    }

    for (auto& th : threads) {
        th.join();
    }

    return result;
}


void Matrix::display() const {
    std::cout << data << std::endl;
}
