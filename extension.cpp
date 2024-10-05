#include <phpcpp.h>
#include "matrix.h"

extern "C" {

PHPCPP_EXPORT void *get_module() {
    static Php::Extension extension("matrix_extension", "1.0");

    Php::Class<Matrix> Matrix("Matrix");

    Matrix.method<&Matrix::__construct>("__construct", {
        Php::ByVal("data", Php::Type::Array, false)
    });

    Matrix.method<&Matrix::add>("add", { Php::ByVal("other", Php::Type::Object) });
    Matrix.method<&Matrix::subtract>("sub", { Php::ByVal("other", Php::Type::Object) });
    Matrix.method<&Matrix::dot>("dot", { Php::ByVal("other", "Matrix") });
    Matrix.method<&Matrix::multiply>("mul", { Php::ByVal("other", Php::Type::Object) });
    Matrix.method<&Matrix::div>("div", {
        Php::ByVal("other", Php::Type::Null, true)
    });

    Matrix.method<&Matrix::exp>("exp", {
        Php::ByVal("other", Php::Type::Numeric, false)
    });

    Matrix.method<&Matrix::inverse>("inverse");
    Matrix.method<&Matrix::determinant>("determinant");
    Matrix.method<&Matrix::log>("log");
    Matrix.method<&Matrix::eigen>("eigen");

    Matrix.method<&Matrix::getItem>("getItem");
    Matrix.method<&Matrix::setItem>("setItem");
    Matrix.method<&Matrix::setData>("setData", { Php::ByVal("data", Php::Type::Array) });

    Matrix.method<&Matrix::getRow>("getRow", {
        Php::ByVal("row", Php::Type::Numeric)
    });

    Matrix.method<&Matrix::setRow>("setRow", {
        Php::ByVal("row", Php::Type::Numeric),
        Php::ByVal("rowData", "Matrix")
    });
    Matrix.method<&Matrix::getCol>("getCol", {
        Php::ByVal("col", Php::Type::Numeric)
    });
    Matrix.method<&Matrix::setCol>("setCol", {
        Php::ByVal("col", Php::Type::Numeric),
        Php::ByVal("colData", "Matrix")
    });

    Matrix.method<&Matrix::argmax>("argmax", {
        Php::ByVal("axis", Php::Type::Numeric, false)
    });

    Matrix.method<&Matrix::clip>("clip", {
        Php::ByVal("min_val", Php::Type::Numeric),
        Php::ByVal("max_val", Php::Type::Numeric)
    });

    Matrix.method<&Matrix::transpose>("transpose");
    Matrix.method<&Matrix::display>("display");
    Matrix.method<&Matrix::getData>("getData");
    Matrix.method<&Matrix::shape>("shape");

    Matrix.method<&Matrix::setThreadScalingFactor>("setThreadScalingFactor", Php::Static, {
        Php::ByVal("factor", Php::Type::Float)
    });

    Matrix.method<&Matrix::random>("random", Php::Static, {
        Php::ByVal("rows", Php::Type::Numeric),
        Php::ByVal("cols", Php::Type::Numeric),
        Php::ByVal("min", Php::Type::Float, false),
        Php::ByVal("max", Php::Type::Float, false),
        Php::ByVal("binomial", Php::Type::Bool, false)
    });


    Matrix.method<&Matrix::round>("round", {
        Php::ByVal("precision", Php::Type::Numeric, false)
    });

    Matrix.method<&Matrix::zeros>("zeros", Php::Static, {
        Php::ByVal("rows", Php::Type::Numeric),
        Php::ByVal("cols", Php::Type::Numeric)
    });

    // Matrix.method<&Matrix::ones_like>("ones_like");

    Matrix.method<&Matrix::getSlice>("getSlice", {
        Php::ByVal("range", Php::Type::Array)
    });

    Matrix.method<&Matrix::zeros_like_static>("zeros_like", Php::Static, {
        Php::ByVal("matrix", "Matrix", false)
    });

    Matrix.method<&Matrix::zeros_like_instance>("zeros_like_instance");

    Matrix.method<&Matrix::ones_like_static>("ones_like", Php::Static, {
        Php::ByVal("matrix", "Matrix", false)
    });

    Matrix.method<&Matrix::ones_like_instance>("ones_like_instance");


    Matrix.method<&Matrix::relu>("relu", {
        Php::ByVal("condition", Php::Type::Numeric),
        Php::ByVal("val", Php::Type::Numeric)
    });

    Matrix.method<&Matrix::offsetGet>("offsetGet", {
        Php::ByVal("key", Php::Type::Numeric)
    });

    Matrix.method<&Matrix::offsetSet>("offsetSet", {
        Php::ByVal("key", Php::Type::Numeric),
        Php::ByVal("value", Php::Type::Array)
    });

    Matrix.method<&Matrix::reshape>("reshape", {
        Php::ByVal("newshape", Php::Type::Array),
        Php::ByVal("order", Php::Type::String, false)
    });

    Matrix.method<&Matrix::sum>("sum", {
        Php::ByVal("a", Php::Type::Numeric, false),
        Php::ByVal("axis", Php::Type::Numeric, false),
        Php::ByVal("initial", Php::Type::Float, false),
        Php::ByVal("where", Php::Type::Array, false)
    });

    Matrix.method<&Matrix::getValuesFromIndices>("getValuesFromIndices", {
        Php::ByVal("indices", Php::Type::Array)
    });

    Matrix.method<&Matrix::pow>("pow", {
        Php::ByVal("exponent", Php::Type::Numeric)
    });
    Matrix.method<&Matrix::abs>("abs");

    Matrix.method<&Matrix::sqrt>("sqrt");

    Matrix.method<&Matrix::eye>("eye", Php::Static, {
        Php::ByVal("N", Php::Type::Numeric),
        Php::ByVal("M", Php::Type::Numeric, false),
        Php::ByVal("k", Php::Type::Numeric, false),
        Php::ByVal("dtype", Php::Type::String, false),
        Php::ByVal("order", Php::Type::String, false)
    });

    Matrix.method<&Matrix::selectRowsByIndices>("selectRowsByIndices", {
        Php::ByVal("indices", Php::Type::Array)
    });

    Matrix.method<&Matrix::mean>("mean", {
        Php::ByVal("axis", Php::Type::Numeric, false),
        Php::ByVal("dtype", Php::Type::String, false),
        Php::ByVal("where", Php::Type::Array, false)
    });

    Matrix.method<&Matrix::std>("std", {
        Php::ByVal("axis", Php::Type::Numeric, false),
        Php::ByVal("ddof", Php::Type::Numeric, false)
    });

    Matrix.method<&Matrix::sign>("sign", { Php::ByVal("where", Php::Type::Array, false) });

    Matrix.method<&Matrix::glorot_uniform>("glorot_uniform", Php::Static, {
        Php::ByVal("fan_in", Php::Type::Numeric),
        Php::ByVal("fan_out", Php::Type::Numeric)
    });

    Matrix.method<&Matrix::min>("min", {
        Php::ByVal("axis", Php::Type::Numeric, true),
        Php::ByVal("initial", Php::Type::Float, true),
        Php::ByVal("where", Php::Type::Array, true)
    });

    Matrix.method<&Matrix::max>("max", {
        Php::ByVal("axis", Php::Type::Numeric, true),
        Php::ByVal("initial", Php::Type::Float, true),
        Php::ByVal("where", Php::Type::Array, true)
    });
    
    extension.add(std::move(Matrix));

    return extension;
}

}