#include "matrixwrapper.h"
#include <iostream>
#include <thread>

void convertPhpArrayToVector(const Php::Value &phpArray, std::vector<std::vector<double>> &outputVector);

void MatrixWrapper::__construct(Php::Parameters &params)
{
    if (params.size() == 1 && params[0].isArray()) {
        Php::Array inputData = params[0];
        int newRows = inputData.size();
        int newCols = newRows > 0 ? Php::count(inputData[0]).numericValue() : 0;

        for (int i = 0; i < newRows; ++i) {
            if (Php::count(inputData[i]).numericValue() != newCols) {
                throw Php::Exception("All rows must have the same number of columns");
            }
        }

        matrix = new Matrix(newRows, newCols);

        for (int i = 0; i < newRows; ++i) {
            Php::Value row = inputData[i];
            if (row.isArray()) {
                Php::Array rowArray = row;
                for (int j = 0; j < newCols; ++j) {
                    (*matrix)(i, j) = rowArray[j];
                }
            } else {
                throw Php::Exception("Each row must be an array");
            }
        }
    } else {
        throw Php::Exception("Invalid parameters for MatrixWrapper constructor");
    }
}

// void MatrixWrapper::__construct(Php::Parameters &params)
// {
//     if (params.size() == 1 && params[0].isArray()) {
//         std::vector<std::vector<double>> inputData;

//         convertPhpArrayToVector(params[0], inputData);

//         int newRows = inputData.size();
//         int newCols = newRows > 0 ? inputData[0].size() : 0;

//         // matrix = new Matrix(inputData);
//         matrix = new Matrix(inputData);
//     } else {
//         throw Php::Exception("Invalid parameters for MatrixWrapper constructor");
//     }
// }



Php::Value MatrixWrapper::add(Php::Parameters &params)
{
    if (params.size() == 1 && params[0].isObject() && params[0].instanceOf("MatrixWrapper")) {
        MatrixWrapper *other = (MatrixWrapper *)params[0].implementation();
        
        Matrix result = matrix->add(*(other->matrix));
        MatrixWrapper *resultWrapper = new MatrixWrapper(std::move(result));
        auto phpObject = Php::Object("MatrixWrapper", resultWrapper);
        
        return phpObject;
    } else if (params.size() == 1) {
        // Scalar addition
        double scalar = params[0];
        Matrix result = matrix->addScalar(scalar);
        return Php::Object("MatrixWrapper", new MatrixWrapper(std::move(result)));
    }

    throw Php::Exception("Invalid parameters for add method");
}


Php::Value MatrixWrapper::div(Php::Parameters& params) {
    if (params[0].isObject()) {
        MatrixWrapper* other = (MatrixWrapper*)params[0].implementation();
        Matrix result = (*matrix) / (*(other->matrix));
        MatrixWrapper* resultWrapper = new MatrixWrapper();
        resultWrapper->matrix = new Matrix(result);
        return Php::Object("MatrixWrapper", resultWrapper);
    } else if (params[0].isNumeric()) {
        double scalar = params[0].numericValue();
        Matrix result = (*matrix) / scalar;
        MatrixWrapper* resultWrapper = new MatrixWrapper();
        resultWrapper->matrix = new Matrix(result);
        return Php::Object("MatrixWrapper", resultWrapper);
    } else {
        throw Php::Exception("Invalid parameter type for division");
    }
}

Php::Value MatrixWrapper::log()
{
    Matrix result = matrix->log();
    return Php::Object("MatrixWrapper", new MatrixWrapper(result));
}

Php::Value MatrixWrapper::exp(Php::Parameters &params)
{
    double base = params.size() > 0 ? params[0].numericValue() : std::exp(1.0); // Default base is e
    Matrix result = matrix->exp(base);
    return Php::Object("MatrixWrapper", new MatrixWrapper(result));
}

Php::Value MatrixWrapper::inverse()
{
    Matrix result = matrix->inverse();
    return Php::Object("MatrixWrapper", new MatrixWrapper(result));
}

Php::Value MatrixWrapper::determinant()
{
    double result = matrix->determinant();
    return Php::Value(result);
}

Php::Value MatrixWrapper::eigen()
{
    auto result = matrix->eigen();
    Php::Array resultArray;
    resultArray[0] = Php::Object("MatrixWrapper", new MatrixWrapper(result.first));
    resultArray[1] = Php::Object("MatrixWrapper", new MatrixWrapper(result.second));
    return resultArray;
}


Php::Value MatrixWrapper::sum(Php::Parameters &params)
{
    int axis = params.size() > 0 ? params[0].numericValue() : -1;
    Matrix result = matrix->sum(axis);
    return Php::Object("MatrixWrapper", new MatrixWrapper(result));
}


Php::Value MatrixWrapper::transpose(){
    Matrix result = matrix->transpose();
    return Php::Object("MatrixWrapper", new MatrixWrapper(result));
}

Php::Value MatrixWrapper::subtract(Php::Parameters &params)
{
    if (params.size() == 1 && params[0].isObject() && params[0].instanceOf("MatrixWrapper")) {
        MatrixWrapper *other = (MatrixWrapper *)params[0].implementation();
        Matrix result = matrix->subtract(*(other->matrix));
        return Php::Object("MatrixWrapper", new MatrixWrapper(result));
    } else if (params.size() == 1) {
        // Scalar subtraction
        double scalar = params[0];
        Matrix result = matrix->subtractScalar(scalar);
        return Php::Object("MatrixWrapper", new MatrixWrapper(result));
    }

    throw Php::Exception("Invalid parameters for subtract method");
}

Php::Value MatrixWrapper::dot(Php::Parameters &params)
{

    if (params[0].isObject()) {
    MatrixWrapper *other = (MatrixWrapper *)params[0].implementation();
    Matrix result = matrix->dot(*(other->matrix));
    return Php::Object("MatrixWrapper", new MatrixWrapper(std::move(result)));
    }else if (params[0].isNumeric()) {
        double scalar = params[0].numericValue();
        Matrix result = (*matrix) * scalar;
        return Php::Object("MatrixWrapper", new MatrixWrapper(std::move(result)));
    } else {
        throw Php::Exception("Invalid parameter type for dot/mul");
    }
}



void MatrixWrapper::setData(Php::Parameters &params)
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
    matrix = new Matrix(newRows, newCols);

    for (int i = 0; i < newRows; ++i) {
        for (int j = 0; j < newCols; ++j) {
            (*matrix)(i, j) = inputData[i][j];
        }
    }
}

Php::Value MatrixWrapper::getData() const
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




Php::Value MatrixWrapper::argmax(Php::Parameters& params) {
    int axis = params.size() > 0 ? params[0].numericValue() : 0;
    std::vector<int> indices = matrix->argmax(axis);

    Php::Array result;
    for (size_t i = 0; i < indices.size(); ++i) {
        result[i] = indices[i];
    }

    return result;
}

Php::Value MatrixWrapper::clip(Php::Parameters &params) {
    double min_val = params[0].numericValue();
    double max_val = params[1].numericValue();
    Matrix result = matrix->clip(min_val, max_val);
    
    // Convert Matrix to Php::Array to return to PHP
    Php::Array phpResult;
    for (int i = 0; i < result.getRows(); ++i) {
        Php::Array row;
        for (int j = 0; j < result.getCols(); ++j) {
            row[j] = result(i, j);
        }
        phpResult[i] = row;
    }
    return phpResult;
}


Php::Value MatrixWrapper::shape(){

    Php::Value array;
    
    array[0] = matrix->getRows();
    
    array[1] = matrix->getCols();


    return array;
}

void MatrixWrapper::display() const
{
    matrix->display();
}

void convertPhpArrayToVector(const Php::Value &phpArray, std::vector<std::vector<double>> &outputVector) {
    int newRows = phpArray.size();
    int newCols = newRows > 0 ? Php::count(phpArray[0]).numericValue() : 0;

    // Resize the output vector
    outputVector.resize(newRows, std::vector<double>(newCols));

    // Initialize the output vector in a single thread
    for (int i = 0; i < newRows; ++i) {
        Php::Value row = phpArray[i];
        if (row.isArray()) {
            Php::Array rowArray = row;
            for (int j = 0; j < newCols; ++j) {
                outputVector[i][j] = rowArray[j];
            }
        } else {
            throw Php::Exception("Each row must be an array");
        }
    }


}