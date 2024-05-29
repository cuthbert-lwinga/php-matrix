#ifndef THREADMANAGER_H
#define THREADMANAGER_H

#include <thread>
#include <vector>
#include <queue>
#include <functional>
#include <future>
#include <mutex>
#include <condition_variable>
#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif

class ThreadManager {
public:
    ThreadManager(size_t numThreads);
    ~ThreadManager();

    template<typename F>
    auto submitTask(F&& task) -> std::future<decltype(task())> {
        auto packagedTask = std::make_shared<std::packaged_task<decltype(task())()>>(std::forward<F>(task));
        std::future<decltype(task())> future = packagedTask->get_future();
        {
            std::unique_lock<std::mutex> lock(queueMutex);
            taskQueue.emplace([packagedTask]() { (*packagedTask)(); });
            ++activeTasks;
        }
        condition.notify_one();
        return future;
    }

    void waitForCompletion();

private:
    void workerThread();
    void gpuWorkerThread();
    void initializeThreads(size_t numThreads);

    std::vector<std::thread> threads;
    std::queue<std::packaged_task<void()>> taskQueue;
    std::mutex queueMutex;
    std::condition_variable condition;
    bool stop;
    size_t activeTasks;
    std::condition_variable completionCondition;

    bool gpuAvailable;

#ifdef USE_CUDA
    cudaStream_t gpuStream;
#endif
};

#endif // THREADMANAGER_H
