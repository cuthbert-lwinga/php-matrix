#include <phpcpp.h>
#include "matrixwrapper.h"

extern "C" {

    PHPCPP_EXPORT void *get_module() {
        static Php::Extension extension("matrix_extension", "1.0");

        Php::Class<MatrixWrapper> matrixWrapper("MatrixWrapper");

        matrixWrapper.method<&MatrixWrapper::__construct>("__construct", {
            Php::ByVal("data", Php::Type::Array,false)
        });

        matrixWrapper.method<&MatrixWrapper::add>("add", { Php::ByVal("other", Php::Type::Object) });
        matrixWrapper.method<&MatrixWrapper::subtract>("sub", { Php::ByVal("other", Php::Type::Object) });
        matrixWrapper.method<&MatrixWrapper::dot>("dot", { Php::ByVal("other", "MatrixWrapper") });
        matrixWrapper.method<&MatrixWrapper::dot>("mul", { Php::ByVal("scalar", Php::Type::Numeric, true) });
        matrixWrapper.method<&MatrixWrapper::div>("div", {
            Php::ByVal("other", Php::Type::Null, true) // Use Php::Type::Null to allow both objects and scalars
        });

        matrixWrapper.method<&MatrixWrapper::log>("exp",{
            Php::ByVal("other", Php::Type::Numeric, true) // Use Php::Type::Null to allow both objects and scalars
        });
        matrixWrapper.method<&MatrixWrapper::log>("sum",{
            Php::ByVal("other", Php::Type::Numeric, true) // Use Php::Type::Null to allow both objects and scalars
        });
        
        matrixWrapper.method<&MatrixWrapper::log>("inverse");
        matrixWrapper.method<&MatrixWrapper::log>("determinant");
        matrixWrapper.method<&MatrixWrapper::log>("log");
        matrixWrapper.method<&MatrixWrapper::log>("eigen");
        matrixWrapper.method<&MatrixWrapper::setData>("setData", { Php::ByVal("data", Php::Type::Array) });

        matrixWrapper.method<&MatrixWrapper::argmax>("argmax", {
            Php::ByVal("axis", Php::Type::Numeric, false)
        });

        matrixWrapper.method<&MatrixWrapper::clip>("clip", {
            Php::ByVal("min_val", Php::Type::Numeric),
            Php::ByVal("max_val", Php::Type::Numeric)
        });

        matrixWrapper.method<&MatrixWrapper::transpose>("transpose");
        matrixWrapper.method<&MatrixWrapper::display>("display");
        matrixWrapper.method<&MatrixWrapper::getData>("getData");
        matrixWrapper.method<&MatrixWrapper::shape>("shape");
        
        matrixWrapper.method<&MatrixWrapper::setThreadScalingFactor>("setThreadScalingFactor", {
            Php::ByVal("factor", Php::Type::Float)
        });


        matrixWrapper.method<&MatrixWrapper::random>("random", {
            Php::ByVal("rows", Php::Type::Numeric),
            Php::ByVal("cols", Php::Type::Numeric),
            Php::ByVal("min", Php::Type::Float, false),
            Php::ByVal("max", Php::Type::Float, false)
        });

        
        // dont use these
        matrixWrapper.method<&MatrixWrapper::offsetGet>("offsetGet", {
            Php::ByVal("key", Php::Type::Numeric)
        });

        matrixWrapper.method<&MatrixWrapper::offsetSet>("offsetSet", {
            Php::ByVal("key", Php::Type::Numeric),
            Php::ByVal("value", Php::Type::Array)
        });



        extension.add(std::move(matrixWrapper));

        return extension;
    }
}
