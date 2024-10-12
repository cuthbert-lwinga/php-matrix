#ifndef MATRIXWRAPPER_H
#define MATRIXWRAPPER_H

#include <Eigen/Dense>
#include <vector>

class MatrixWrapper {

public:
    Eigen::MatrixXd data;
    static double threadScalingFactor;
    int threads = 8;

    MatrixWrapper(int rows, int cols, double value = 0.0);
    MatrixWrapper(const std::vector<std::vector<double>> &inputData);
    MatrixWrapper(const Eigen::MatrixXd& other);

    static MatrixWrapper glorot_uniform(int fan_in, int fan_out);

    // Static method to set the thread scaling factor
    static void setThreadScalingFactor(double factor) {
        threadScalingFactor = factor;
    }
    // Getters

    double getItem(int row, int col) const;
    void setItem(int row, int col, double value);

    int getRows() const { return data.rows(); }
    int getCols() const { return data.cols(); }

    void setRow(int row, const MatrixWrapper& rowData);
    MatrixWrapper getCol(int col) const;
    void setCol(int col, const MatrixWrapper& colData);

    MatrixWrapper getRow(int row) const;

    // Arithimetic Operations (addition, subtraction, multiplication, division)

    MatrixWrapper add(const MatrixWrapper &other) const;
    MatrixWrapper subtract(const MatrixWrapper &other) const;
    MatrixWrapper dot(const MatrixWrapper &other) const;
    MatrixWrapper multiply(const MatrixWrapper &other) const;
    MatrixWrapper multiplyElementwise(const MatrixWrapper &other) const;
    MatrixWrapper addScalar(const double other) const;
    MatrixWrapper subtractScalar(const double other) const;
    MatrixWrapper log() const;
    MatrixWrapper sqrt() const; 
    MatrixWrapper exp(double base = std::exp(1.0)) const;
    MatrixWrapper sum(int axis = -1) const;
    MatrixWrapper std(int axis = -1, double ddof = 0) const;
    MatrixWrapper max(int axis = -1, double initial = std::numeric_limits<double>::lowest(), const Eigen::MatrixXd& where = Eigen::MatrixXd::Constant(0, 0, 1.0)) const;
    MatrixWrapper min(int axis = -1, double initial = std::numeric_limits<double>::max(), const Eigen::MatrixXd& where = Eigen::MatrixXd::Constant(0, 0, 1.0)) const;

    MatrixWrapper getSlice(const std::vector<int>& range) const;
    static MatrixWrapper zeros_like(const MatrixWrapper& other);
    MatrixWrapper zeros_like() const;

    static MatrixWrapper ones_like(const MatrixWrapper& other);
    MatrixWrapper ones_like() const;

    MatrixWrapper copy() const;

    MatrixWrapper round(int precision);


    MatrixWrapper inverse() const;
    double determinant() const;
    std::pair<MatrixWrapper, MatrixWrapper> eigen() const; // Eigenvalues and Eigenvectors

    // Other functions    
    std::vector<int> argmax(int axis = 0) const;
    MatrixWrapper clip(double min_val, double max_val) const;
    MatrixWrapper transpose() const; 
    static MatrixWrapper random(int rows, int cols, double min = 0.0, double max = 1.0);
    MatrixWrapper randomBinomial(int rows, int cols);
    MatrixWrapper relu(double condition, double val) const;
    MatrixWrapper reshape(const std::vector<int>& newshape, const char order = 'C') const;    
    MatrixWrapper abs() const;
    MatrixWrapper sum(int axis = -1, double initial = 0.0, const Eigen::MatrixXd& where = Eigen::MatrixXd::Constant(0, 0, 0.0)) const;

    MatrixWrapper pow(const MatrixWrapper& exponent) const;
    MatrixWrapper pow(double exponent) const;
    MatrixWrapper getValuesFromIndices(const std::vector<int>& indices) const;
    static MatrixWrapper eye(int N, int M = -1, int k = 0, const std::string& dtype = "float", const std::string& order = "C");
    MatrixWrapper selectRowsByIndices(const std::vector<int>& indices) const;
    MatrixWrapper mean(int axis = -1, const std::string& dtype = "float64", const Eigen::MatrixXd& where = Eigen::MatrixXd::Constant(0, 0, 0.0)) const;
    MatrixWrapper sign(const Eigen::MatrixXd& where = Eigen::MatrixXd::Constant(0, 0, 1.0)) const;
    MatrixWrapper selectByIndices(const std::vector<int>& rowIndices, const std::vector<int>& colIndices) const;
    MatrixWrapper slice(int start, int length, int axis = 1) const;


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
