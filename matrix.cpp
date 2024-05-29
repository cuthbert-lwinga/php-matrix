#include "matrix.h"
#include <iostream>
#include <thread>
#include <vector>
#include <map> 
#include <algorithm>
#include <cmath>
// ThreadManager MatrixWrapper::threadManager(4);
int calculateThreadsBasedOnMatrixSize__(int rows, int cols, double scalingFactor, int maxThreads);

void Matrix::__construct(Php::Parameters &params)
{

   if (params.size() == 1 && params[0].isArray()) {
        Php::Value inputData = params[0];
        std::vector<std::vector<double>> data;

        for (auto &row : inputData) {
            std::vector<double> rowData;
            for (auto &element : row.second) {
                rowData.push_back(element.second.floatValue());
            }
            data.push_back(rowData);
        }

        int rows = data.size();
        int cols = rows > 0 ? data[0].size() : 0;

        // Initialize the Matrix object with the input data using parallelization
        matrix = new MatrixWrapper(rows, cols);


        // Initialize the matrix in parallel using Eigen's internal parallelization
        #pragma omp parallel for
        for (int i = 0; i < rows; ++i) {
            matrix->data.row(i) = Eigen::Map<const Eigen::RowVectorXd>(&data[i][0], cols);
        }
    } else {
        
        matrix = new MatrixWrapper(1, 1);


        //throw Php::Exception("Invalid parameters for MatrixWrapper constructor");
    }



}


Php::Value Matrix::add(Php::Parameters &params)
{
    if (params.size() == 1 && params[0].isObject() && params[0].instanceOf("Matrix")) {
        Matrix *other = (Matrix *)params[0].implementation();
        MatrixWrapper result = matrix->add(*(other->matrix));
        auto phpObject = Php::Object("Matrix", new Matrix(std::move(result)));
        return phpObject;

    } else if (params.size() == 1) {
        // Scalar addition
        double scalar = params[0];
        MatrixWrapper result = matrix->addScalar(scalar);
        return Php::Object("Matrix", new Matrix(std::move(result)));
    }

    throw Php::Exception("Invalid parameters for add method");
}


Php::Value Matrix::div(Php::Parameters &params)
{
    if (params.size() == 1 && params[0].isObject() && params[0].instanceOf("Matrix")) {
        
        Matrix *other = (Matrix *)params[0].implementation();
        MatrixWrapper result = (*matrix) / (*(other->matrix));
        auto phpObject = Php::Object("Matrix", new Matrix(std::move(result)));
        return phpObject;

    } else if (params.size() == 1) {
        // Scalar addition
        double scalar = params[0];
        MatrixWrapper result = (*matrix) / scalar;

        return Php::Object("Matrix", new Matrix(std::move(result)));
    }

    throw Php::Exception("Invalid parameters for add method");
}


Php::Value Matrix::log()
{
    MatrixWrapper result = matrix->log();
    return Php::Object("Matrix", new Matrix(result));
}

Php::Value Matrix::exp(Php::Parameters &params)
{
    double base = params.size() > 0 ? params[0].numericValue() : std::exp(1.0); // Default base is e
    MatrixWrapper result = matrix->exp(base);
    return Php::Object("Matrix", new Matrix(result));
}

Php::Value Matrix::inverse()
{
    MatrixWrapper result = matrix->inverse();
    return Php::Object("Matrix", new Matrix(result));
}

Php::Value Matrix::determinant()
{
    double result = matrix->determinant();
    return Php::Value(result);
}

Php::Value Matrix::eigen()
{
    auto result = matrix->eigen();
    Php::Array resultArray;
    resultArray[0] = Php::Object("Matrix", new Matrix(result.first));
    resultArray[1] = Php::Object("Matrix", new Matrix(result.second));
    return resultArray;
}


Php::Value Matrix::sum(Php::Parameters &params)
{
    int axis = params.size() > 0 ? params[0].numericValue() : -1;
    MatrixWrapper result = matrix->sum(axis);
    return Php::Object("Matrix", new Matrix(result));
}


Php::Value Matrix::transpose(){
    MatrixWrapper result = matrix->transpose();
    return Php::Object("Matrix", new Matrix(result));
}

Php::Value Matrix::subtract(Php::Parameters &params)
{
    if (params.size() == 1 && params[0].isObject() && params[0].instanceOf("Matrix")) {
        Matrix *other = (Matrix *)params[0].implementation();
        MatrixWrapper result = matrix->subtract(*(other->matrix));
        return Php::Object("Matrix", new Matrix(result));
    } else if (params.size() == 1) {
        // Scalar subtraction
        double scalar = params[0];
        MatrixWrapper result = matrix->subtractScalar(scalar);
        return Php::Object("Matrix", new Matrix(result));
    }

    throw Php::Exception("Invalid parameters for subtract method");
}

Php::Value Matrix::dot(Php::Parameters &params)
{

    if (params[0].isObject()) {
        Matrix *other = (Matrix *)params[0].implementation();
        MatrixWrapper result = matrix->dot(*(other->matrix));
        return Php::Object("Matrix", new Matrix(std::move(result)));
    }else if (params[0].isNumeric()) {
        double scalar = params[0].numericValue();
        MatrixWrapper result = (*matrix) * scalar;
        return Php::Object("Matrix", new Matrix(std::move(result)));
    } else {
        throw Php::Exception("Invalid parameter type for dot/mul");
    }
}



void Matrix::setData(Php::Parameters &params)
{
    Php::Array inputData = params[0];
    int newRows = inputData.size();
    int newCols = newRows > 0 ? ((Php::Array)inputData[0]).size() : 0;

    for (int i = 0; i < newRows; ++i) {
        if (((Php::Array)inputData[i]).size() != newCols) {
            throw Php::Exception("All rows must have the same number of columns");
        }
    }

    delete matrix;
    matrix = new MatrixWrapper(newRows, newCols);

    for (int i = 0; i < newRows; ++i) {
        for (int j = 0; j < newCols; ++j) {
            (*matrix)(i, j) = inputData[i][j];
        }
    }
}

Php::Value Matrix::getData() const
{
    Php::Array outputData;
    for (int i = 0; i < matrix->getRows(); ++i) {
        Php::Array row;
        for (int j = 0; j < matrix->getCols(); ++j) {
            row[j] = (*matrix)(i, j);
        }
        outputData[i] = row;
    }
    return outputData;
}




Php::Value Matrix::argmax(Php::Parameters& params) {
    int axis = params.size() > 0 ? params[0].numericValue() : 0;
    std::vector<int> indices = matrix->argmax(axis);

    Php::Array result;
    for (size_t i = 0; i < indices.size(); ++i) {
        result[i] = indices[i];
    }

    return result;
}

Php::Value Matrix::clip(Php::Parameters &params) {
    double min_val = params[0].numericValue();
    double max_val = params[1].numericValue();
    MatrixWrapper result = matrix->clip(min_val, max_val);
    
    auto phpObject = Php::Object("Matrix", new Matrix(std::move(result)));
    return phpObject;
}


Php::Value Matrix::shape(){

    Php::Value array;
    
    array[0] = matrix->getRows();
    
    array[1] = matrix->getCols();


    return array;
}


Php::Value Matrix::random(Php::Parameters &params) {
    int rows = params[0].numericValue();
    int cols = params[1].numericValue();
    double min = params.size() > 2 ? params[2].floatValue() : 0.0;
    double max = params.size() > 3 ? params[3].floatValue() : 1.0;

    MatrixWrapper result = matrix->random(rows, cols, min, max);

    matrix = new MatrixWrapper(std::move(result));

    auto phpObject = Php::Object("Matrix", new Matrix(std::move(result)));
    return phpObject;
}


void Matrix::display() const
{
    matrix->display();
}


Php::Value Matrix::offsetGet(Php::Parameters &params)
{
    int row = params[0].numericValue();
    if (params.size() == 2) {
        int col = params[1].numericValue();
        return (*matrix)(row, col);
    }
    // Return a row as an array
    std::vector<double> row_data;
    for (int col = 0; col < matrix->getCols(); ++col) {
        row_data.push_back((*matrix)(row, col));
    }
    return row_data;
}

void Matrix::offsetSet(Php::Parameters &params) {
    if (params.size() != 2 || !params[0].isNumeric() || !params[1].isArray()) {
        throw Php::Exception("Invalid argument(s) for offsetSet");
    }

    int index = params[0].numericValue();

    if (index >= matrix->getRows() || index < 0) {
        throw Php::Exception("Index out of bounds");
    }

    Php::Array phpArray = params[1];

    if (phpArray.size() != matrix->getCols()) {
        throw Php::Exception("Array size does not match matrix column count");
    }

    Eigen::VectorXd vec(matrix->getCols());
    for (int i = 0; i < phpArray.size(); ++i) {
        vec(i) = phpArray[i];
    }

    matrix->data.row(index) = vec;
}

// Helper functions

// Function to calculate a safe number of threads for context switching based on matrix size
int calculateThreadsBasedOnMatrixSize__(int rows, int cols, double scalingFactor = 1.0, int maxThreads = 2500) {
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

