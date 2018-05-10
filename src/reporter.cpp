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

bool Reporter::isNumeric(std::string s)
{
    if (!s.empty() && s[0] != '-') s = "0" + s; // prepend 0
    std::string garbage;
    std::stringstream ss(s);
    ss >> *(std::unique_ptr<double>(new double)) >> garbage;
    return garbage.empty();
}

void Reporter::report(const std::string& input, bool wrap)
{
    // split by new line
    std::vector<std::string> paragraph;
    split(paragraph, input.c_str(), '\n');
    std::string error_prefix("Error:");
    std::string warning_prefix("Warning:");
    std::string error_color_start = "\033[1;31m";
    std::string warning_color_start = "\033[1;33m";
    std::string color_end = "\033[0m";
    std::string list_prefix;
    for (auto&& message : paragraph) {
        std::vector<std::string> line;
        // red for error, yellow for warnings
        bool is_error = !message.compare(0, error_prefix.size(), error_prefix);
        bool is_warning =
            !message.compare(0, warning_prefix.size(), warning_prefix);
        bool is_list = false;
        if (!is_error && !is_warning) {
            std::vector<std::string> check_list;
            split(check_list, message.c_str(), ')');
            if (check_list.size() > 1) {
                if (isNumeric(check_list.front())) {
                    is_list = true;
                    list_prefix = check_list.front() + ")";
                }
            }
        }
        split(line, message.c_str(), ' ');
        if (!wrap) {
            if (is_error) std::cerr << error_color_start;
            if (is_warning) std::cerr << warning_color_start;
            std::cerr << message << '\n';
            if (is_error || is_warning) std::cerr << color_end;
            m_log_file << message << '\n';
        }
        else
        {
            size_t cur_length = 0;
            if (is_error) std::cerr << error_color_start;
            if (is_warning) std::cerr << warning_color_start;
            for (auto&& word : line) {
                if (word.length() + cur_length >= m_width) {
                    // word is too long, so display on the next line
                    cur_length = word.length() + 1;
                    std::cerr << '\n';
                    m_log_file << '\n';
                    if (is_error) {
                        std::cerr
                            << std::string(error_prefix.length() + 1, ' ');
                        m_log_file
                            << std::string(error_prefix.length() + 1, ' ');
                    }
                    else if (is_warning)
                    {
                        std::cerr
                            << std::string(warning_prefix.length() + 1, ' ');
                        m_log_file
                            << std::string(warning_prefix.length() + 1, ' ');
                    }
                    else if (is_list)
                    {
                        std::cerr << std::string(list_prefix.length() + 1, ' ');
                        m_log_file
                            << std::string(list_prefix.length() + 1, ' ');
                    }
                }
                else
                {
                    cur_length += word.length() + 1;
                }
                std::cerr << word << " ";
                m_log_file << word << " ";
            }
            if (is_error || is_warning) std::cerr << color_end;
            std::cerr << '\n';
            m_log_file << '\n';
        }
    }
    std::cerr << '\n' << std::flush;
    m_log_file << '\n' << std::flush;
}
