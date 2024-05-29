#include <iostream>
#include <chrono>

int main() {
    const int numIterations = 1000000000;
    auto start = std::chrono::high_resolution_clock::now();

    // Loop over 1,000,000 times
    for (int i = 0; i < numIterations; ++i) {
        // Simple operation
        int x = i * 2;
        //std::cout<<"=> "<<i<<std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Time taken to loop over " << numIterations << " times: " << elapsed.count() << " seconds" << std::endl;
    return 0;
}
