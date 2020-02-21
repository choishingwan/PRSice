// This file is part of PRSice-2, copyright (C) 2016-2019
// Shing Wan Choi, Paul F. O’Reilly
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef THREAD_QUEUE_H
#define THREAD_QUEUE_H

#include <condition_variable>
#include <mutex>
#include <queue>
#include <thread>

template <typename T>
class Thread_Queue
{
public:
    bool pop(T& item)
    {
        bool completed = false;
        std::unique_lock<std::mutex> mlock(m_mutex);
        m_cond_not_empty.wait(
            mlock, [this] { return (m_storage_queue.size() || m_completed); });
        completed = m_completed;
        if (!completed)
        {
            item = std::move(m_storage_queue.front());
            m_storage_queue.pop();
        }
        m_num_processing--;
        mlock.unlock();
        m_cond_not_full.notify_one();
        return completed;
    }

    bool pop(T& item, size_t num_thread)
    {
        std::unique_lock<std::mutex> mlock(m_mutex);
        m_cond_not_empty.wait(mlock, [this, num_thread] {
            return (m_storage_queue.size() || (m_num_completed == num_thread));
        });
        if (m_num_completed != num_thread)
        {
            item = std::move(m_storage_queue.front());
            m_storage_queue.pop();
        }
        m_num_processing--;
        mlock.unlock();
        m_cond_not_full.notify_one();
        return (m_num_completed == num_thread);
    }
    void push(const T& item, size_t max_process)
    {
        std::unique_lock<std::mutex> mlock(m_mutex);
        // stop producer from producing extra data when
        // we have not finish enough jobs
        while (max_process <= m_num_processing) { m_cond_not_full.wait(mlock); }
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
        while (max_process <= m_num_processing) { m_cond_not_full.wait(mlock); }
        m_storage_queue.push(std::move(item));
        m_num_processing++;
        mlock.unlock();
        m_cond_not_empty.notify_one();
    }
    void emplace(T&& item, size_t max_process)
    {
        std::unique_lock<std::mutex> mlock(m_mutex);
        while (max_process <= m_num_processing) { m_cond_not_full.wait(mlock); }
        m_storage_queue.emplace(std::forward<T>(item));
        // m_storage_queue.push(item);
        m_num_processing++;
        mlock.unlock();
        m_cond_not_empty.notify_one();
    }
    void emplace(T&& item)
    {
        std::unique_lock<std::mutex> mlock(m_mutex);
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
        ++m_num_completed;
        mlock.unlock();
        m_cond_not_empty.notify_all();
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
    size_t m_num_completed = 0;
    bool m_completed = false;
};


#endif
