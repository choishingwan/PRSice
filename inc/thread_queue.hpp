
#ifndef THREAD_QUEUE_H
#define THREAD_QUEUE_H

#include <queue>
#ifdef _WIN32
#include <mingw.condition_variable.h>
#include <mingw.mutex.h>
#include <mingw.thread.h>
#else
#include <condition_variable>
#include <mutex>
#include <thread>
#endif

template <typename T>
class Thread_Queue
{
public:
    void pop(T& item)
    {
        std::unique_lock<std::mutex> mlock(m_mutex);
        while (m_storage_queue.empty()) {
            m_cond_not_empty.wait(mlock);
        }
        item = std::move(m_storage_queue.front());
        m_num_processing--;
        m_storage_queue.pop();
        mlock.unlock();
        m_cond_not_full.notify_one();
    }

    void push(const T& item, size_t max_process)
    {
        std::unique_lock<std::mutex> mlock(m_mutex);
        // stop producer from producing extra data when
        // we have not finish enough jobs
        while (max_process <= m_num_processing) {
            m_cond_not_full.wait(mlock);
        }
        m_storage_queue.push(item);
        m_num_processing++;
        mlock.unlock();
        m_cond_not_empty.notify_one();
    }
    void push(T&& item, size_t max_process)
    {
        std::unique_lock<std::mutex> mlock(m_mutex);
        // stop producer from producing extra data when
        // we have not finish enough jobs
        while (max_process <= m_num_processing) {
            m_cond_not_full.wait(mlock);
        }
        m_storage_queue.push(item);
        m_num_processing++;
        mlock.unlock();
        m_cond_not_empty.notify_one();
    }
    void emplace(T&& item, size_t max_process)
    {
        std::unique_lock<std::mutex> mlock(m_mutex);
        while (max_process <= m_num_processing) {
            m_cond_not_full.wait(mlock);
        }
        m_storage_queue.emplace(std::forward<T>(item));
        // m_storage_queue.push(item);
        m_num_processing++;
        mlock.unlock();
        m_cond_not_empty.notify_one();
    }
    void completed()
    {
        std::unique_lock<std::mutex> mlock(m_mutex);
        m_completed = true;
        mlock.unlock();
        m_cond_not_empty.notify_one();
    }
    bool has_completed()
    {
        bool completed = false;
        std::unique_lock<std::mutex> mlock(m_mutex);
        completed = m_completed;
        mlock.unlock();
        m_cond_not_full.notify_one();
        return completed;
    }
    size_t num_processing() const { return m_num_processing; }
    Thread_Queue() = default;
    Thread_Queue(const Thread_Queue&) = delete;            // disable copying
    Thread_Queue& operator=(const Thread_Queue&) = delete; // disable assignment
private:
    std::queue<T> m_storage_queue;
    std::mutex m_mutex;
    std::condition_variable m_cond_not_empty;
    std::condition_variable m_cond_not_full;
    size_t m_num_processing = 0;
    bool m_completed = false;
};


#endif
