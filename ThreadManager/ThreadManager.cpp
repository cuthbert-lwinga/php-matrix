#include "ThreadManager.h"
#include <iostream>
#include <system_error>

ThreadManager::ThreadManager(size_t numThreads)
    : stop(false), activeTasks(0)
{
#ifdef USE_CUDA
    // Check for GPU availability
    int deviceCount;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);
    gpuAvailable = (err == cudaSuccess && deviceCount > 0);
    if (gpuAvailable) {
        cudaStreamCreate(&gpuStream);
    }
#endif

    initializeThreads(numThreads);
}

void ThreadManager::initializeThreads(size_t numThreads) {
    bool initialized = false;

    // Try the user-specified number of threads first
    size_t threadsToInitialize = numThreads;

    while (!initialized && threadsToInitialize > 0) {
        try {
            // Clear any existing threads
            threads.clear();

            // Initialize worker threads
            for (size_t i = 0; i < threadsToInitialize; ++i) {
                threads.emplace_back(&ThreadManager::workerThread, this);
            }

#ifdef USE_CUDA
            if (gpuAvailable) {
                threads.emplace_back(&ThreadManager::gpuWorkerThread, this);
            }
#endif

            initialized = true;
        } catch (const std::system_error& e) {
            std::cerr << "Failed to initialize ThreadManager with " << threadsToInitialize << " threads: " << e.what() << std::endl;
            threadsToInitialize /= 2;
        }
    }

    // If initialization failed with the specified number of threads, fall back to hardware concurrency
    if (!initialized) {
        threadsToInitialize = std::thread::hardware_concurrency();
        std::cerr << "Falling back to hardware concurrency: " << threadsToInitialize << " threads." << std::endl;
        while (!initialized && threadsToInitialize > 0) {
            try {
                // Clear any existing threads
                threads.clear();

                // Initialize worker threads
                for (size_t i = 0; i < threadsToInitialize; ++i) {
                    threads.emplace_back(&ThreadManager::workerThread, this);
                }

#ifdef USE_CUDA
                if (gpuAvailable) {
                    threads.emplace_back(&ThreadManager::gpuWorkerThread, this);
                }
#endif

                initialized = true;
            } catch (const std::system_error& e) {
                std::cerr << "Failed to initialize ThreadManager with " << threadsToInitialize << " threads: " << e.what() << std::endl;
                threadsToInitialize /= 2;
            }
        }
    }

    // If initialization still failed, fall back to a single-threaded setup
    if (!initialized) {
        threadsToInitialize = 1;
        std::cerr << "Falling back to single-threaded setup." << std::endl;
        try {
            // Clear any existing threads
            threads.clear();

            // Initialize worker thread
            threads.emplace_back(&ThreadManager::workerThread, this);

#ifdef USE_CUDA
            if (gpuAvailable) {
                threads.emplace_back(&ThreadManager::gpuWorkerThread, this);
            }
#endif

        } catch (const std::system_error& e) {
            throw std::runtime_error("Could not initialize ThreadManager even with a single thread.");
        }
    }
}

ThreadManager::~ThreadManager()
{
    {
        std::unique_lock<std::mutex> lock(queueMutex);
        stop = true;
    }
    condition.notify_all();
    for (std::thread& thread : threads) {
        if (thread.joinable()) {
            thread.join();
        }
    }

#ifdef USE_CUDA
    if (gpuAvailable) {
        cudaStreamDestroy(gpuStream);
    }
#endif
}

void ThreadManager::waitForCompletion()
{
    std::unique_lock<std::mutex> lock(queueMutex);
    completionCondition.wait(lock, [this] { return activeTasks == 0; });
}

void ThreadManager::workerThread()
{
    while (true) {
        std::packaged_task<void()> task;
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            condition.wait(lock, [this] { return !taskQueue.empty() || stop; });
            if (stop && taskQueue.empty()) return;
            task = std::move(taskQueue.front());
            taskQueue.pop();
        }
        task();
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            --activeTasks;
            if (activeTasks == 0) {
                completionCondition.notify_all();
            }
        }
    }
}

#ifdef USE_CUDA
void ThreadManager::gpuWorkerThread()
{
    while (true) {
        std::packaged_task<void()> task;
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            condition.wait(lock, [this] { return !taskQueue.empty() || stop; });
            if (stop && taskQueue.empty()) return;
            task = std::move(taskQueue.front());
            taskQueue.pop();
        }
        // Execute the task on GPU
        task();
        cudaStreamSynchronize(gpuStream);
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            --activeTasks;
            if (activeTasks == 0) {
                completionCondition.notify_all();
            }
        }
    }
}
#endif
