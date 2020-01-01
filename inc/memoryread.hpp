#ifndef MEMORYREAD_HPP
#define MEMORYREAD_HPP

#include "misc.hpp"
#include <fstream>
#include <mio.hpp>
#include <stdexcept>
#include <string>


class FileRead
{
public:
    FileRead() {}
    void read(const std::string& file, const std::streampos& byte_pos,
              const std::streampos read_size, char* result)
    {
        if (file != m_file_name) { new_file(file, byte_pos); }

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
            throw std::runtime_error("Error: Cannot read file: " + m_file_name);
        }
        m_offset = read_size + byte_pos;
    }

private:
    std::ifstream m_input;
    std::string m_file_name;
    std::streampos m_offset;
    void new_file(const std::string& file, const std::streampos byte_pos)
    {
        m_file_name = file;
        m_offset = byte_pos;
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
};

#endif // MEMORYREAD_HPP
