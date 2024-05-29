CXX = g++
CXXFLAGS = -std=c++14 -pthread

# Uncomment the following lines if using CUDA
# CUDA_PATH = /usr/local/cuda
# CXXFLAGS += -DUSE_CUDA -I$(CUDA_PATH)/include
# LDFLAGS = -L$(CUDA_PATH)/lib64 -lcudart

TARGET = thread_manager_test
SOURCES = main.cpp ThreadManager.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJECTS)

.PHONY: all clean