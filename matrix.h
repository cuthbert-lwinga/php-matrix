#ifndef MATRIX_H
#define MATRIX_H

#include <Eigen/Dense>
#include <vector>

class Matrix {

public:
    Eigen::MatrixXd data;
    static double threadScalingFactor;
    int threads = 10;

    Matrix(int rows, int cols, double value = 0.0);
    Matrix(const std::vector<std::vector<double>> &inputData);
    Matrix(const Eigen::MatrixXd& other);


    // Static method to set the thread scaling factor
    static void setThreadScalingFactor(double factor) {
        threadScalingFactor = factor;
    }
    // Getters

    int getRows() const { return data.rows(); }
    int getCols() const { return data.cols(); }

    // Arithimetic Operations (addition, subtraction, multiplication, division)

    Matrix add(const Matrix &other) const;
    Matrix subtract(const Matrix &other) const;
    Matrix dot(const Matrix &other) const;
    Matrix addScalar(const double other) const;
    Matrix subtractScalar(const double other) const;
    Matrix log() const;
    Matrix sqrt() const; 
    Matrix exp(double base = std::exp(1.0)) const;
    Matrix sum(int axis = -1) const;
    Matrix inverse() const;
    double determinant() const;
    std::pair<Matrix, Matrix> eigen() const; // Eigenvalues and Eigenvectors

    // Other functions    
    std::vector<int> argmax(int axis = 0) const;
    Matrix clip(double min_val, double max_val) const;
    Matrix transpose() const; 
    Matrix random(int rows, int cols, double min = 0.0, double max = 1.0);


    // Operator overloads for +, -, *
    double& operator()(int row, int col) { return data(row, col); }
    const double& operator()(int row, int col) const { return data(row, col); }

    Matrix operator+(const Matrix &other) const {
        return add(other);
    }

    Matrix operator-(const Matrix &other) const ;
    Matrix operator-(double scalar) const;

    Matrix operator*(const Matrix &other) const;
    Matrix operator*(double scalar) const;
    // Inside the Matrix class definition
    Matrix operator/(const Matrix& other) const;
    Matrix operator/(double scalar) const;


    void display() const;
};


#endif // MATRIX_H
