EXTENSION_NAME=matrix_extension
PHP_CONFIG=$(shell which php-config)
PHPCPP_DIR=/usr/local/include
PHPCPP_LIB=/usr/local/lib
PHP_EXTENSION_DIR=$(shell $(PHP_CONFIG) --extension-dir)
EIGEN_DIR=./eigen

CXX=g++
CXXFLAGS=-O3 -Wall -c -std=c++14 -fpic `$(PHP_CONFIG) --includes` -I$(PHPCPP_DIR) -I$(EIGEN_DIR)
LDFLAGS=-shared -L$(PHPCPP_LIB) -lphpcpp

SRC=matrix.cpp matrixwrapper.cpp extension.cpp
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

install: $(EXTENSION_NAME).so
	cp $(EXTENSION_NAME).so $(PHP_EXTENSION_DIR)

clean:
	rm -f $(OBJ) $(EXTENSION_NAME).so

.PHONY: all install clean
