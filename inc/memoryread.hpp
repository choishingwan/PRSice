#ifndef MEMORYREAD_HPP
#define MEMORYREAD_HPP

#include "misc.hpp"
#include <fstream>
#include <mio.hpp>
#include <stdexcept>
#include <string>

// TODO: Might want to do some precomputation of the ideal offsets such that we
// can minimize the number of time we do the mapping
class MemoryRead
{
public:
    MemoryRead() {}
    void read(const std::string& file, const unsigned long long& byte_pos,
              const unsigned long long read_size, char* result)
    {
        if (file != m_file_name) { new_file(file, byte_pos); }
        if (m_use_mmap)
        {
            // first check if we need to remap the file
            // last condition is for bgen, which require read of different size
            // will still be suboptimal for bgen unless we know exactly how many
            // byte each data span as a whole
            if (byte_pos >= (m_offset + m_block_size) || byte_pos < m_offset
                || byte_pos - m_offset + read_size > m_block_size)
            {
                // + byte_pos to account for possible useless bytes
                m_offset = byte_pos;
                std::error_code error;
                m_memory_map.map(m_file_name, m_offset, m_block_size, error);
                if (error)
                {
                    throw std::runtime_error("Error: Failed to map file: "
                                             + m_file_name);
                }
            }
            if (byte_pos - m_offset + read_size > m_memory_map.mapped_length())
            {
                // As we have re-mapped by this point, the only possible reason
                // for over-run is read_size > file size
                throw std::runtime_error("Error: File read out of bound");
            }
            char* read = result;
            for (unsigned long long i = 0; i < read_size; ++i)
            {
                *read = m_memory_map[byte_pos - m_offset + i];
                ++read;
            }
        }
        else
        {
            assert(m_input.is_open());
            if (byte_pos != m_offset)
            {
                if (!m_input.seekg(byte_pos, std::ios_base::beg))
                {
                    throw std::runtime_error("Error: Cannot seek within file: "
                                             + m_file_name);
                }
            }
            if (!m_input.read(result, read_size))
            {
                throw std::runtime_error("Error: Cannot read file: "
                                         + m_file_name);
            }
            m_offset = read_size + byte_pos;
        }
    }
    void init_memory_map(const unsigned long long mem,
                         const unsigned long long& data_size)
    {
        bool allow_mmap = calculate_block_size(mem, data_size);
        if (m_use_mmap) { m_use_mmap = allow_mmap; }
        if (!allow_mmap)
        {
            std::cerr << "Warning: Not enough memory for file mapping to be "
                         "worth it, will "
                         "fall back to traditional file read"
                      << std::endl;
        }
        m_mem_calculated = true;
    }
    bool mem_calculated() const { return m_mem_calculated; }
    void no_mmap() { m_use_mmap = false; }
    void use_mmap() { m_use_mmap = true; }

private:
    std::ifstream m_input;
    mio::mmap_source m_memory_map;
    std::string m_file_name;
    unsigned long long m_offset;
    unsigned long long m_block_size;
    bool m_use_mmap = false;
    bool m_mem_calculated = false;
    // return true if we think mmap is useful
    bool calculate_block_size(const unsigned long long& mem,
                              const unsigned long long& data_size)
    {
        unsigned long long remain_mem = misc::remain_memory();
        if (mem > remain_mem)
        {
            std::cerr << "Warning: Not enough memory left. Only " << remain_mem
                      << " byte remain, but requested " << mem
                      << " byte. Will fall back to tradition file read"
                      << std::endl;
            return false;
        }
        unsigned long long num_block = (mem / data_size);
        m_block_size = num_block * data_size;
        return (num_block > 1);
    }
    void new_file(const std::string& file, const unsigned long long byte_pos)
    {
        m_file_name = file;
        m_offset = byte_pos;
        if (m_use_mmap)
        {
            m_offset = byte_pos;
            std::error_code error;
            auto&& file = mio::detail::open_file(m_file_name,
                                                 mio::access_mode::read, error);
            if (mio::detail::query_file_size(file, error) < byte_pos)
            {
                // when we have enough memory to read the whole file, do that
                m_offset = 0;
            }

            m_memory_map.map(m_file_name, m_offset, m_block_size, error);
            if (error)
            {
                throw std::runtime_error("Error: Cannot map file: "
                                         + m_file_name);
            }
        }
        else
        {
            if (m_input.is_open()) { m_input.close(); }
            m_input.clear();
            m_input.open(m_file_name.c_str(), std::ios::binary);
            if (byte_pos != 0)
            {
                if (!m_input.seekg(byte_pos, std::ios_base::beg))
                {
                    throw std::runtime_error("Error: Cannot seek within file: "
                                             + m_file_name);
                }
            }
        }
    }
};

#endif // MEMORYREAD_HPP
