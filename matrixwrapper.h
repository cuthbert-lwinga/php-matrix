#ifndef MATRIXWRAPPER_H
#define MATRIXWRAPPER_H

#include <Eigen/Dense>
#include <vector>

class MatrixWrapper {

public:
    Eigen::MatrixXd data;
    static double threadScalingFactor;
    int threads = 10;

    MatrixWrapper(int rows, int cols, double value = 0.0);
    MatrixWrapper(const std::vector<std::vector<double>> &inputData);
    MatrixWrapper(const Eigen::MatrixXd& other);


    // Static method to set the thread scaling factor
    static void setThreadScalingFactor(double factor) {
        threadScalingFactor = factor;
    }
    // Getters

    int getRows() const { return data.rows(); }
    int getCols() const { return data.cols(); }

    // Arithimetic Operations (addition, subtraction, multiplication, division)

    MatrixWrapper add(const MatrixWrapper &other) const;
    MatrixWrapper subtract(const MatrixWrapper &other) const;
    MatrixWrapper dot(const MatrixWrapper &other) const;
    MatrixWrapper addScalar(const double other) const;
    MatrixWrapper subtractScalar(const double other) const;
    MatrixWrapper log() const;
    MatrixWrapper sqrt() const; 
    MatrixWrapper exp(double base = std::exp(1.0)) const;
    MatrixWrapper sum(int axis = -1) const;
    MatrixWrapper inverse() const;
    double determinant() const;
    std::pair<MatrixWrapper, MatrixWrapper> eigen() const; // Eigenvalues and Eigenvectors

    // Other functions    
    std::vector<int> argmax(int axis = 0) const;
    MatrixWrapper clip(double min_val, double max_val) const;
    MatrixWrapper transpose() const; 
    MatrixWrapper random(int rows, int cols, double min = 0.0, double max = 1.0);


    // Operator overloads for +, -, *
    double& operator()(int row, int col) { return data(row, col); }
    const double& operator()(int row, int col) const { return data(row, col); }

    MatrixWrapper operator+(const MatrixWrapper &other) const {
        return add(other);
    }

    MatrixWrapper operator-(const MatrixWrapper &other) const ;
    MatrixWrapper operator-(double scalar) const;

    MatrixWrapper operator*(const MatrixWrapper &other) const;
    MatrixWrapper operator*(double scalar) const;
    // Inside the Matrix class definition
    MatrixWrapper operator/(const MatrixWrapper& other) const;
    MatrixWrapper operator/(double scalar) const;


    void display() const;
};


#endif // MATRIX_H
