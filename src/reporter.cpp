/*
 * Reporter.cpp
 *
 *  Created on: 1 Nov 2017
 *      Author: shingwan
 */

#include <reporter.hpp>

void split(std::vector<std::string>& result, const char* str, char c = ' ')
{
    do
    {
        const char* begin = str;

        while (*str != c && *str) str++;

        result.push_back(std::string(begin, str));
    } while (0 != *str++);
}


void Reporter::report(const std::string& input, bool wrap)
{
    // split by new line
    std::vector<std::string> paragraph;
    split(paragraph, input.c_str(), '\n');
    for (auto&& message : paragraph) {
        std::vector<std::string> line;
        split(line, message.c_str(), ' ');
        if (!wrap) {
            std::cerr << message << '\n';
            m_log_file << message << '\n';
        }
        else
        {
            size_t cur_length = 0;
            for (auto&& word : line) {
                if (word.length() + cur_length >= m_width) {
                    cur_length = word.length() + 1;
                    std::cerr << '\n';
                    m_log_file << '\n';
                }
                else
                {
                    cur_length += word.length() + 1;
                }
                std::cerr << word << " ";
                m_log_file << word << " ";
            }
            std::cerr << '\n';
            m_log_file << '\n';
        }
    }
    std::cerr << '\n' << std::flush;
    m_log_file << '\n' << std::flush;
}
