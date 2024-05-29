#ifndef MATRIXWRAPPER_H
#define MATRIXWRAPPER_H

#include <phpcpp.h>
#include "matrix.h"
#include "ThreadManager.h"

class MatrixWrapper : public Php::Base
{
private:
    Matrix* matrix;

public:
    MatrixWrapper() {}
    MatrixWrapper(const Matrix &matrix) : matrix(new Matrix(matrix)) {}
    virtual ~MatrixWrapper() { delete matrix; }
    // static ThreadManager threadManager;

    void __construct(Php::Parameters &params);
    void setData(Php::Parameters &params);
    
    // Arithimetic Operations (addition, subtraction, multiplication, division)
    Php::Value add(Php::Parameters &params);
    Php::Value subtract(Php::Parameters &params);
    Php::Value div(Php::Parameters& params);
    Php::Value dot(Php::Parameters &params);
    Php::Value log();
    Php::Value exp(Php::Parameters &params); // Update to accept parameters
    Php::Value sum(Php::Parameters &params);
    Php::Value inverse();
    Php::Value determinant();
    Php::Value eigen();


    // Other functions
    Php::Value argmax(Php::Parameters& params);
    Php::Value transpose();
    Php::Value shape();
    Php::Value clip(Php::Parameters &params);
    Php::Value random(Php::Parameters &params);
    
    Php::Value offsetGet(Php::Parameters &params);
    void offsetSet(Php::Parameters &params);

    Php::Value getData() const;
    void display() const;


    // Method to set the thread scaling factor
    static void setThreadScalingFactor(Php::Parameters &params) {
        double factor = params[0].floatValue();
        Matrix::setThreadScalingFactor(factor);
    }

};

#endif // MATRIXWRAPPER_H
