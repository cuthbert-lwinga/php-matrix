EXTENSION_NAME=matrix_extension

PHP_CONFIG=$(shell which php-config)
<<<<<<< HEAD

PHPCPP_DIR=$(shell $(PHP_CONFIG) --includes | sed 's/-I/-I\/usr\/include\/phpcpp /g')

PHPCPP_LIB=/usr/local/lib

PHP_EXTENSION_DIR=$(shell $(PHP_CONFIG) --extension-dir)

EIGEN_DIR=./eigen
=======
PHPCPP_DIR=/usr/include/phpcpp
PHPCPP_LIB=/usr/lib
PHP_EXTENSION_DIR=$(shell $(PHP_CONFIG) --extension-dir)
EIGEN_DIR=./Eigen
PHP_CONFIG_DIR=/etc/php/8.1/cli/conf.d/
PHPINIFILE=30-phpmatrix.ini
>>>>>>> 262795e (local changhes)

THREAD_MANAGER_DIR=./ThreadManager

CXX=g++

CXXFLAGS=-O3 -Wall -c -g -std=c++14 -fpic `$(PHP_CONFIG) --includes` -I$(PHPCPP_DIR) -I$(EIGEN_DIR) -I$(THREAD_MANAGER_DIR)

LDFLAGS=-shared -L$(PHPCPP_LIB) -lphpcpp -lpthread

SRC=matrix.cpp matrixwrapper.cpp extension.cpp $(THREAD_MANAGER_DIR)/ThreadManager.cpp

OBJ=$(SRC:.cpp=.o)

all: $(EXTENSION_NAME).so

$(EXTENSION_NAME).so: $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $^

matrix.o: matrix.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

matrixwrapper.o: matrixwrapper.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

extension.o: extension.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

ThreadManager.o: $(THREAD_MANAGER_DIR)/ThreadManager.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

install: $(EXTENSION_NAME).so
	cp $(EXTENSION_NAME).so $(PHP_EXTENSION_DIR)

clean:
	rm -f $(OBJ) $(EXTENSION_NAME).so

.PHONY: all install clean