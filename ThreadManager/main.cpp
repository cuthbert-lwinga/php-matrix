#include "ThreadManager.h"
#include <iostream>
#include <vector>
#include <mutex>
#include <chrono>

std::mutex printMutex;

void testTask(int chunkId, int rowStart, int rowEnd, int colStart, int colEnd, std::vector<std::vector<int>>& chunkMatrix, int threadId) {
    {
        std::lock_guard<std::mutex> lock(printMutex);
        std::cout << "Task for chunk " << chunkId << " by thread " << threadId << " is starting." << std::endl;
    }

    for (int row = rowStart; row < rowEnd; ++row) {
        for (int col = colStart; col < colEnd; ++col) {
            chunkMatrix[row - rowStart][col - colStart] = threadId; // Assign the thread id to the matrix element
        }
    }

    std::this_thread::sleep_for(std::chrono::seconds(1));

    {
        std::lock_guard<std::mutex> lock(printMutex);
        std::cout << "Task for chunk " << chunkId << " by thread " << threadId << " is finished." << std::endl;
    }
}

int main() {
    int rows, cols, numThreads;
    std::cout << "Enter number of rows: ";
    std::cin >> rows;
    std::cout << "Enter number of columns: ";
    std::cin >> cols;
    std::cout << "Enter number of threads: ";
    std::cin >> numThreads;

    ThreadManager threadManager(numThreads);

    // Define chunk size
    int chunkSize = 10000; // Adjust this value based on available memory

    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();

    // Process the matrix in chunks
    int chunkId = 0;
    for (int rowStart = 0; rowStart < rows; rowStart += chunkSize) {
        for (int colStart = 0; colStart < cols; colStart += chunkSize) {
            int rowEnd = std::min(rowStart + chunkSize, rows);
            int colEnd = std::min(colStart + chunkSize, cols);

            // Create the chunk matrix
            std::vector<std::vector<int>> chunkMatrix(rowEnd - rowStart, std::vector<int>(colEnd - colStart, 0));

            // Submit task to process the chunk
            threadManager.submitTask([chunkId, rowStart, rowEnd, colStart, colEnd, &chunkMatrix, threadId=chunkId % numThreads] {
                testTask(chunkId, rowStart, rowEnd, colStart, colEnd, chunkMatrix, threadId);
            });

            // Wait for the chunk to be processed before moving to the next chunk
            threadManager.waitForCompletion();

            // Process the chunk as needed (e.g., store or output the chunk)
            // Example: std::cout << "Chunk " << chunkId << " processed.\n";
            chunkId++;
        }
    }

    // Stop the timer
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "All tasks are complete." << std::endl;
    std::cout << "Time taken: " << elapsed.count() << " seconds" << std::endl;

    return 0;
}
