#ifndef MATRIX_H
#define MATRIX_H

#include <phpcpp.h>
#include "matrixwrapper.h"
#include "ThreadManager/ThreadManager.h"

class Matrix : public Php::Base
{
private:
    MatrixWrapper matrix;

public:
    // Default constructor: create a 1x1 matrix with value 0
    Matrix() : matrix(1, 1, 0.0) {}
    
    // Copy constructor
    Matrix(const MatrixWrapper &temp) : matrix(temp) {}
    
    // Constructor with dimensions
    Matrix(int rows, int cols) : matrix(rows, cols) {}
    
    virtual ~Matrix() {}  
    // static ThreadManager threadManager;

    void __construct(Php::Parameters &params);
    void setData(Php::Parameters &params);
    
    Php::Value getItem(Php::Parameters &params);
    void setItem(Php::Parameters &params);
    
    Php::Value getRow(Php::Parameters &params);
    void setRow(Php::Parameters &params);
    Php::Value getCol(Php::Parameters &params);
    void setCol(Php::Parameters &params);
    
    // Arithimetic Operations (addition, subtraction, multiplication, division)
    Php::Value add(Php::Parameters &params);
    Php::Value subtract(Php::Parameters &params);
    Php::Value div(Php::Parameters& params);
    Php::Value dot(Php::Parameters &params);
    Php::Value multiply(Php::Parameters &params);
    Php::Value log();
    Php::Value exp(Php::Parameters &params); // Update to accept parameters
    Php::Value sum(Php::Parameters &params);
    Php::Value inverse();
    Php::Value determinant();
    Php::Value eigen();
    Php::Value selectByIndices(Php::Parameters &params);


    // Other functions
    Php::Value max(Php::Parameters &params);
    Php::Value min(Php::Parameters &params);
    Php::Value argmax(Php::Parameters& params);
    Php::Value transpose();
    Php::Value shape();
    Php::Value clip(Php::Parameters &params);
    static Php::Value random(Php::Parameters &params);
    static Php::Value zeros(Php::Parameters &params);
    //Php::Value ones_like(Php::Parameters &params);

    Php::Value getSlice(Php::Parameters& params);
    static Php::Value zeros_like_static(Php::Parameters& params);
    Php::Value zeros_like_instance();

    static Php::Value ones_like_static(Php::Parameters& params);
    Php::Value ones_like_instance();

    Php::Value relu(Php::Parameters &params);
    Php::Value reshape(Php::Parameters& params);
    Php::Value abs();
    Php::Value offsetGet(Php::Parameters &params);
    Php::Value sqrt();
    Php::Value pow(Php::Parameters &params);
    Php::Value getValuesFromIndices(Php::Parameters &params);
    static Php::Value eye(Php::Parameters &params);
    Php::Value selectRowsByIndices(Php::Parameters &params);
    Php::Value mean(Php::Parameters &params);
    Php::Value sign(Php::Parameters &params);
    Php::Value std(Php::Parameters &params);
    Php::Value round(Php::Parameters &params);
    static Php::Value glorot_uniform(Php::Parameters &params);
    static Php::Value oneHotEncoded(Php::Parameters &params);
    Php::Value slice(Php::Parameters &params);



    Php::Value copy();

    void offsetSet(Php::Parameters &params);


    Php::Value getData() const;
    void display() const;


    // Method to set the thread scaling factor
    static void setThreadScalingFactor(Php::Parameters &params) {
        double factor = params[0].floatValue();
        MatrixWrapper::setThreadScalingFactor(factor);
    }

};

#endif // MATRIXWRAPPER_H
