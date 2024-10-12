#include "matrix.h"
#include <iostream>
#include <thread>
#include <vector>
#include <map> 
#include <algorithm>
#include <cmath>
#include <random>

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
        
        matrix = new MatrixWrapper(1, 1,0.0);


        //throw Php::Exception("Invalid parameters for MatrixWrapper constructor");
    }



}

Php::Value Matrix::getItem(Php::Parameters &params)
{
    if (!matrix) {
        throw Php::Exception("Matrix is not initialized");
    }
    int row = params[0];
    int col = params[1];
    try {
        return matrix->getItem(row, col);
    } catch (const std::out_of_range& e) {
        throw Php::Exception(e.what());
    }
}

void Matrix::setItem(Php::Parameters &params)
{
    if (!matrix) {
        throw Php::Exception("Matrix is not initialized");
    }
    int row = params[0];
    int col = params[1];
    double value = params[2];
    try {
        matrix->setItem(row, col, value);
    } catch (const std::out_of_range& e) {
        throw Php::Exception(e.what());
    }
}

Php::Value Matrix::getRow(Php::Parameters &params)
{
    int row = params[0];
    try {
        MatrixWrapper rowMatrix = matrix->getRow(row);
        Matrix* newMatrix = new Matrix(rowMatrix);
        return Php::Object("Matrix", newMatrix);
    }
    catch (const std::out_of_range& e) {
        throw Php::Exception(e.what());
    }
}

void Matrix::setRow(Php::Parameters &params)
{
    int row = params[0];
    Php::Value rowData = params[1];
    if (!rowData.instanceOf("Matrix")) {
        throw Php::Exception("Row data must be a Matrix object");
    }
    Matrix *rowMatrix = (Matrix *)rowData.implementation();
    try {
        matrix->setRow(row, *(rowMatrix->matrix));
    }
    catch (const std::exception& e) {
        throw Php::Exception(e.what());
    }
}

Php::Value Matrix::getCol(Php::Parameters &params)
{
    int col = params[0];
    try {
        MatrixWrapper colMatrix = matrix->getCol(col);
        Matrix* newMatrix = new Matrix(colMatrix);
        return Php::Object("Matrix", newMatrix);
    }
    catch (const std::exception& e) {
        throw Php::Exception(e.what());
    }
}

void Matrix::setCol(Php::Parameters &params)
{
    int col = params[0];
    Php::Value colData = params[1];
    if (!colData.instanceOf("Matrix")) {
        throw Php::Exception("Column data must be a Matrix object");
    }
    Matrix *colMatrix = (Matrix *)colData.implementation();
    try {
        matrix->setCol(col, *(colMatrix->matrix));
    }
    catch (const std::exception& e) {
        throw Php::Exception(e.what());
    }
}

Php::Value Matrix::add(Php::Parameters &params) {
    if (params.size() == 1 && params[0].isObject() && params[0].instanceOf("Matrix")) {
        Matrix *other = (Matrix *)params[0].implementation();
        MatrixWrapper result = matrix->add(*(other->matrix));
        return Php::Object("Matrix", new Matrix(result));
    } else if (params.size() == 1) {
        double scalar = params[0].floatValue();
        MatrixWrapper result = matrix->addScalar(scalar);
        return Php::Object("Matrix", new Matrix(result));
    }

    throw Php::Exception("Invalid parameters for add method");
}


Php::Value Matrix::div(Php::Parameters &params) {
    if (params.empty()) {
        throw Php::Exception("Division requires an argument");
    }

    if (params[0].isNumeric()) {
        // Scalar division
        double scalar = params[0].numericValue();
        return Php::Object("Matrix", new Matrix(*matrix / scalar));
    } else if (params[0].isObject() && params[0].instanceOf("Matrix")) {
        // Matrix division
        Matrix *other = (Matrix *)params[0].implementation();
        
        try {
            MatrixWrapper result = *matrix / *(other->matrix);
            return Php::Object("Matrix", new Matrix(result));
        } catch (const std::exception& e) {
            throw Php::Exception(e.what());
        }
    } else {
        throw Php::Exception("Invalid argument type for division");
    }
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
        double scalar = params[0].floatValue();
        MatrixWrapper result = matrix->subtractScalar(scalar);
        return Php::Object("Matrix", new Matrix(result));
    }

    throw Php::Exception("Invalid parameters for subtract method");
}


Php::Value Matrix::dot(Php::Parameters &params)
{
    if (params.size() == 1) {
        if (params[0].isObject() && params[0].instanceOf("Matrix")) {
            Matrix *other = (Matrix *)params[0].implementation();
            MatrixWrapper result = matrix->dot(*(other->matrix));
            return Php::Object("Matrix", new Matrix(result));
        } else {
            double scalar = params[0].numericValue();
            MatrixWrapper result = (*matrix) * scalar;
            return Php::Object("Matrix", new Matrix(result));
        } 
    } else {
        throw Php::Exception("Invalid number of parameters for dot/mul");
    }
}


Php::Value Matrix::multiply(Php::Parameters &params) {
    if (params.size() == 1) {
        if (params[0].isObject() && params[0].instanceOf("Matrix")) {
            Matrix *other = (Matrix *)params[0].implementation();
            MatrixWrapper result(matrix->getRows(), matrix->getCols());
            if (matrix->getCols() == other->matrix->getCols() && matrix->getRows() == other->matrix->getRows()) {
                result = matrix->multiplyElementwise(*(other->matrix));
            } else if (other->matrix->getCols() == 1 && matrix->getRows() == other->matrix->getRows()) {
                result = matrix->multiplyElementwise(*(other->matrix));
            } else {
                throw Php::Exception("Matrix dimensions must agree for element-wise multiplication or allow for broadcasting.");
            }
            return Php::Object("Matrix", new Matrix(result));
        } else {
            double scalar = params[0].floatValue();
            MatrixWrapper result = (*matrix) * scalar;
            return Php::Object("Matrix", new Matrix(result));
        }
    } else {
        throw Php::Exception("Invalid number of parameters for multiply");
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

// In matrix.cpp
Php::Value Matrix::min(Php::Parameters &params)
{
    int axis = -1;
    double initial = std::numeric_limits<double>::max();
    Eigen::MatrixXd where;

    if (params.size() > 0) axis = params[0].numericValue();
    if (params.size() > 1) initial = params[1].numericValue();
    if (params.size() > 2 && params[2].isArray()) {
        Php::Value phpWhere = params[2];
        where = Eigen::MatrixXd::Constant(matrix->getRows(), matrix->getCols(), 1.0);
        for (int i = 0; i < matrix->getRows(); ++i) {
            for (int j = 0; j < matrix->getCols(); ++j) {
                where(i, j) = phpWhere[i][j] ? 1.0 : 0.0;
            }
        }
    }

    MatrixWrapper result = matrix->min(axis, initial, where);
    return Php::Object("Matrix", new Matrix(result));
}

Php::Value Matrix::clip(Php::Parameters &params) {
    double min_val = params[0];
    double max_val = params[1];
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

Php::Value Matrix::selectByIndices(Php::Parameters &params) {
    if (params.size() != 2 || !params[0].isArray() || !params[1].isArray()) {
        throw Php::Exception("selectByIndices expects two array parameters");
    }

    std::vector<int> rowIndices;
    std::vector<int> colIndices;

    for (const auto& value : params[0]) {
        rowIndices.push_back(value.second.numericValue());
    }

    for (const auto& value : params[1]) {
        colIndices.push_back(value.second.numericValue());
    }

    MatrixWrapper result = matrix->selectByIndices(rowIndices, colIndices);
    return Php::Object("Matrix", new Matrix(result));
}

Php::Value Matrix::random(Php::Parameters &params) {
    int rows = params[0].numericValue();
    int cols = params[1].numericValue();
    double min = params.size() > 2 ? params[2].floatValue() : 0.0;
    double max = params.size() > 3 ? params[3].floatValue() : 1.0;
    bool binomial = params.size() > 4 ? params[4].boolValue() : false;

    // MatrixWrapper result(1, 1);
    // if (binomial) {
    //     result = matrix->randomBinomial(rows, cols);
    // } else {
    //     result = matrix->random(rows, cols, min, max);
    // }

    MatrixWrapper matrix = MatrixWrapper::random(rows, cols, min, max);

    // delete matrix;
    // matrix = new MatrixWrapper(result);

    return Php::Object("Matrix", new Matrix(matrix));
}

Php::Value Matrix::round(Php::Parameters &params) {
    int precision = 2;
    if (params.size() == 1 && params[0].isNumeric()) {
       precision = params[0].numericValue();
    }


    MatrixWrapper result = matrix->round(precision);
    return Php::Object("Matrix", new Matrix(result));
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


Php::Value Matrix::zeros(Php::Parameters &params)
{
    int rows = params[0];
    int cols = params[1];
    MatrixWrapper result(rows, cols,0.0);
    auto phpObject = Php::Object("Matrix", new Matrix(std::move(result)));
    return phpObject;
    
    }



Php::Value Matrix::getSlice(Php::Parameters& params) {
    if (params.size() != 1 || !params[0].isArray()) {
        throw Php::Exception("getSlice expects an array parameter");
    }
    std::vector<int> range;
    for (const auto& value : params[0]) {
        range.push_back(value.second.numericValue());
    }
    MatrixWrapper result = matrix->getSlice(range);
    return Php::Object("Matrix", new Matrix(result));
}

Php::Value Matrix::zeros_like_static(Php::Parameters& params) {
    if (params.size() > 0 && params[0].isObject()) {
        Matrix* other = (Matrix*)params[0].implementation();
        MatrixWrapper result = MatrixWrapper::zeros_like(*(other->matrix));
        return Php::Object("Matrix", new Matrix(result));
    }
    throw Php::Exception("zeros_like expects a Matrix object parameter");
}

Php::Value Matrix::ones_like_instance() {
    MatrixWrapper result = matrix->ones_like();
    return Php::Object("Matrix", new Matrix(result));
}

Php::Value Matrix::ones_like_static(Php::Parameters& params) {
    if (params.size() > 0 && params[0].isObject()) {
        Matrix* other = (Matrix*)params[0].implementation();
        MatrixWrapper result = MatrixWrapper::ones_like(*(other->matrix));
        return Php::Object("Matrix", new Matrix(result));
    }
    throw Php::Exception("zeros_like expects a Matrix object parameter");
}

Php::Value Matrix::zeros_like_instance() {
    MatrixWrapper result = matrix->zeros_like();
    return Php::Object("Matrix", new Matrix(result));
}

Php::Value Matrix::relu(Php::Parameters &params) {
    if (params.size() != 2) {
        throw Php::Exception("Invalid number of parameters for relu");
    }

    double condition = params[0];
    double val = params[1];

    MatrixWrapper result = matrix->relu(condition, val);
    return Php::Object("Matrix", new Matrix(std::move(result)));
}

// matrix.cpp
Php::Value Matrix::reshape(Php::Parameters& params) {
    if (params.size() != 1 && params.size() != 2) {
        throw Php::Exception("Invalid number of parameters for reshape");
    }

    Php::Value newshape = params[0];
    if (!newshape.isArray()) {
        throw Php::Exception("First parameter must be an array representing the new shape");
    }

    char order = params.size() == 2 ? params[1].stringValue()[0] : 'C';
    if (order != 'C' && order != 'F') {
        throw Php::Exception("Invalid order parameter. Must be 'C' or 'F'");
    }

    std::vector<int> shape_vec;
    for (auto& value : newshape) {
        shape_vec.push_back(value.second.numericValue());
    }

    MatrixWrapper result = matrix->reshape(shape_vec, order);
    return Php::Object("Matrix", new Matrix(result));
}


Php::Value Matrix::abs() {
    MatrixWrapper result = matrix->abs();
    return Php::Object("Matrix", new Matrix(result));
}


Php::Value Matrix::sum(Php::Parameters &params) {
    int axis = params.size() > 0 ? params[0].numericValue() : -1;
    double initial = params.size() > 1 ? params[1].floatValue() : 0.0;
    Eigen::MatrixXd where = Eigen::MatrixXd::Constant(matrix->getRows(), matrix->getCols(), 0.0);

    if (params.size() > 2 && params[2].isArray()) {
        Php::Array phpArray = params[2];
        for (int i = 0; i < matrix->getRows(); ++i) {
            for (int j = 0; j < matrix->getCols(); ++j) {
                where(i, j) = phpArray[i][j] ? 1.0 : 0.0;
            }
        }
    }

    MatrixWrapper result = matrix->sum(axis, initial, where);
    return Php::Object("Matrix", new Matrix(result));
}


Php::Value Matrix::sqrt() {
    MatrixWrapper result = matrix->sqrt();
    return Php::Object("Matrix", new Matrix(result));
}

// Php::Array MatrixWrapper::getData() const {
//     Php::Array outputData;
//     for (int i = 0; i < data.rows(); ++i) {
//         Php::Array row;
//         for (int j = 0; j < data.cols(); ++j) {
//             row[j] = data(i, j);
//         }
//         outputData[i] = row;
//     }
//     return outputData;
// }

Php::Value Matrix::pow(Php::Parameters &params) {
    if (params.size() == 1 && params[0].isObject() && params[0].instanceOf("Matrix")) {
        Matrix* exponent = (Matrix*)params[0].implementation();
        MatrixWrapper result = matrix->pow(*exponent->matrix);
        return Php::Object("Matrix", new Matrix(result));
    } else if (params.size() == 1 && params[0].isNumeric()) {
        double exponent = params[0].numericValue();
        MatrixWrapper result = matrix->pow(exponent);
        return Php::Object("Matrix", new Matrix(result));
    } else {
        throw Php::Exception("Invalid parameters for pow method");
    }
}


Php::Value Matrix::getValuesFromIndices(Php::Parameters &params) {
    if (params.size() != 1 || !params[0].isArray()) {
        throw Php::Exception("Invalid parameter for getValuesFromIndices");
    }

    std::vector<int> indices;
    Php::Array phpArray = params[0];
    for (auto &value : phpArray) {
        indices.push_back(value.second);
    }

    MatrixWrapper result = matrix->getValuesFromIndices(indices);
    return Php::Object("Matrix", new Matrix(result));
}

Php::Value Matrix::eye(Php::Parameters &params) {
    int N = params[0].numericValue();
    int M = params.size() > 1 ? params[1].numericValue() : -1;
    int k = params.size() > 2 ? params[2].numericValue() : 0;
    std::string dtype = params.size() > 3 ? params[3].stringValue() : "float";
    std::string order = params.size() > 4 ? params[4].stringValue() : "C";

    MatrixWrapper result = MatrixWrapper::eye(N, M, k, dtype, order);
    return Php::Object("Matrix", new Matrix(result));
}


Php::Value Matrix::selectRowsByIndices(Php::Parameters &params) {
    if (params.size() != 1 || !params[0].isArray()) {
        throw Php::Exception("Invalid parameter for selectRowsByIndices");
    }

    std::vector<int> indices;
    Php::Array phpArray = params[0];
    for (auto &value : phpArray) {
        indices.push_back(value.second);
    }

    MatrixWrapper result = matrix->selectRowsByIndices(indices);
    return Php::Object("Matrix", new Matrix(result));
}


Php::Value Matrix::mean(Php::Parameters &params) {
    int axis = params.size() > 0 ? params[0].numericValue() : -1;
    std::string dtype = params.size() > 1 ? params[1].stringValue() : "float64";
    Eigen::MatrixXd where = Eigen::MatrixXd::Constant(matrix->getRows(), matrix->getCols(), 0.0);

    if (params.size() > 2 && params[2].isArray()) {
        Php::Array phpArray = params[2];
        for (int i = 0; i < matrix->getRows(); ++i) {
            for (int j = 0; j < matrix->getCols(); ++j) {
                where(i, j) = phpArray[i][j] ? 1.0 : 0.0;
            }
        }
    }

    MatrixWrapper result = matrix->mean(axis, dtype, where);
    return Php::Object("Matrix", new Matrix(result));
}

Php::Value Matrix::sign(Php::Parameters &params) {
    Eigen::MatrixXd where = Eigen::MatrixXd::Constant(matrix->getRows(), matrix->getCols(), 1.0);

    if (params.size() > 0 && params[0].isArray()) {
        Php::Array phpArray = params[0];
        for (int i = 0; i < matrix->getRows(); ++i) {
            for (int j = 0; j < matrix->getCols(); ++j) {
                where(i, j) = phpArray[i][j] ? 1.0 : 0.0;
            }
        }
    }

    MatrixWrapper result = matrix->sign(where);
    return Php::Object("Matrix", new Matrix(result));
}


Php::Value Matrix::std(Php::Parameters &params) {
    int axis = params.size() > 0 ? params[0].numericValue() : -1; // Default to flatten array
    double ddof = params.size() > 1 ? params[1].numericValue() : 0; // Default ddof is 0
    MatrixWrapper result = matrix->std(axis, ddof);
    return Php::Object("Matrix", new Matrix(result));
}


Php::Value Matrix::copy()
{
    if (!matrix) {
        throw Php::Exception("Matrix is not initialized");
    }
    Matrix* newMatrix = new Matrix(matrix->copy());
    return Php::Object("Matrix", newMatrix);
}


Php::Value Matrix::oneHotEncoded(Php::Parameters &params) {
    if (params.size() < 2) {
        throw Php::Exception("oneHotEncoded requires at least two parameters: numRows and indices");
    }

    int numRows = params[0];
    Php::Value indices = params[1];
    
    // Manually find the maximum index
    int maxIndex = -1;
    for (const auto &index : indices) {
        int currentIndex = index.second.numericValue();
        if (currentIndex > maxIndex) {
            maxIndex = currentIndex;
        }
    }

    int numCols = params.size() > 2 ? params[2].numericValue() : maxIndex + 1;

    MatrixWrapper result(numRows, numCols, 0.0);

    int row = 0;
    for (const auto &index : indices) {
        int col = index.second.numericValue();
        if (col >= 0 && col < numCols) {
            result.setItem(row, col, 1.0);
        }
        row++;
    }

    return Php::Object("Matrix", new Matrix(result));
}

Php::Value Matrix::max(Php::Parameters &params)
{
    int axis = -1;
    double initial = std::numeric_limits<double>::lowest();
    Eigen::MatrixXd where;

    if (params.size() > 0) axis = params[0].numericValue();
    if (params.size() > 1) initial = params[1].numericValue();
    if (params.size() > 2 && params[2].isArray()) {
        Php::Value phpWhere = params[2];
        where = Eigen::MatrixXd::Constant(matrix->getRows(), matrix->getCols(), 1.0);
        for (int i = 0; i < matrix->getRows(); ++i) {
            for (int j = 0; j < matrix->getCols(); ++j) {
                where(i, j) = phpWhere[i][j] ? 1.0 : 0.0;
            }
        }
    }

    MatrixWrapper result = matrix->max(axis, initial, where);
    return Php::Object("Matrix", new Matrix(result));
}


Php::Value Matrix::glorot_uniform(Php::Parameters &params)
{
    if (params.size() != 2) {
        throw Php::Exception("Glorot uniform initialization requires 2 parameters: fan_in and fan_out");
    }

    int fan_in = params[0].numericValue();
    int fan_out = params[1].numericValue();

    MatrixWrapper result = MatrixWrapper::glorot_uniform(fan_in, fan_out);
    return Php::Object("Matrix", new Matrix(result));
}

Php::Value Matrix::slice(Php::Parameters &params) {
    if (params.size() < 2) {
        throw Php::Exception("Slice requires at least start and length parameters");
    }

    int start = params[0].numericValue();
    int length = params[1].numericValue();
    int axis = params.size() > 2 ? params[2].numericValue() : 0;

    try {
        MatrixWrapper result = matrix->slice(start, length, axis);
        return Php::Object("Matrix", new Matrix(result));
    } catch (const std::invalid_argument& e) {
        throw Php::Exception(e.what());
    }
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

