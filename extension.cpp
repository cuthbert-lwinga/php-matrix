#include <phpcpp.h>
#include "matrix.h"

extern "C" {

    PHPCPP_EXPORT void *get_module() {
        static Php::Extension extension("matrix_extension", "1.0");

        Php::Class<Matrix> Matrix("Matrix");

        Matrix.method<&Matrix::__construct>("__construct", {
            Php::ByVal("data", Php::Type::Array,false)
        });

        Matrix.method<&Matrix::add>("add", { Php::ByVal("other", Php::Type::Object) });
        Matrix.method<&Matrix::subtract>("sub", { Php::ByVal("other", Php::Type::Object) });
        Matrix.method<&Matrix::dot>("dot", { Php::ByVal("other", "Matrix") });
        Matrix.method<&Matrix::dot>("mul", { Php::ByVal("scalar", Php::Type::Numeric, true) });
        Matrix.method<&Matrix::div>("div", {
            Php::ByVal("other", Php::Type::Null, true) // Use Php::Type::Null to allow both objects and scalars
        });

        Matrix.method<&Matrix::log>("exp",{
            Php::ByVal("other", Php::Type::Numeric, true) // Use Php::Type::Null to allow both objects and scalars
        });
        Matrix.method<&Matrix::log>("sum",{
            Php::ByVal("other", Php::Type::Numeric, true) // Use Php::Type::Null to allow both objects and scalars
        });
        
        Matrix.method<&Matrix::log>("inverse");
        Matrix.method<&Matrix::log>("determinant");
        Matrix.method<&Matrix::log>("log");
        Matrix.method<&Matrix::log>("eigen");
        Matrix.method<&Matrix::setData>("setData", { Php::ByVal("data", Php::Type::Array) });

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
        
        Matrix.method<&Matrix::setThreadScalingFactor>("setThreadScalingFactor", {
            Php::ByVal("factor", Php::Type::Float)
        });


        Matrix.method<&Matrix::random>("random", {
            Php::ByVal("rows", Php::Type::Numeric),
            Php::ByVal("cols", Php::Type::Numeric),
            Php::ByVal("min", Php::Type::Float, false),
            Php::ByVal("max", Php::Type::Float, false)
        });

        
        // dont use these
        Matrix.method<&Matrix::offsetGet>("offsetGet", {
            Php::ByVal("key", Php::Type::Numeric)
        });

        Matrix.method<&Matrix::offsetSet>("offsetSet", {
            Php::ByVal("key", Php::Type::Numeric),
            Php::ByVal("value", Php::Type::Array)
        });



        extension.add(std::move(Matrix));

        return extension;
    }
}
