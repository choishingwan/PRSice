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

#ifndef REPORTER_HPP_
#define REPORTER_HPP_

#include "misc.hpp"
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

class Reporter
{
public:
    Reporter() {}
    Reporter(bool test) : m_unit_test(test) {}
    Reporter(const std::string& log_name, size_t width = 60, bool test = false)
        : m_width(width), m_unit_test(test)
    {
        if (!test)
        {
            m_log_file.open(log_name.c_str());
            if (!m_log_file.is_open())
            {
                throw std::runtime_error("Error: " + log_name
                                         + " cannot be open");
            }
        }
    }
    void initialize(const std::string& log_name, size_t width = 60)
    {
        m_width = width;
        if (m_log_file.is_open()) m_log_file.close();
        m_log_file.open(log_name.c_str());
        if (!m_log_file.is_open())
        {
            throw std::runtime_error("Error: " + log_name + " cannot be open");
        }
    }
    virtual ~Reporter();
    void report(const std::string& input, bool wrap = true);
    void simple_report(const std::string& input);

private:
    std::ofstream m_log_file;
    const std::string m_error_prefix = "Error:";
    const std::string m_warning_prefix = "Warning:";
    const size_t m_error_prefix_size = 6;
    const size_t m_warning_prefix_size = 8;
    size_t m_width = 60;
    bool m_unit_test = false;
#if defined(WIN32) || defined(_WIN32) \
    || defined(__WIN32) && !defined(__CYGWIN__)
    // some windows might support colored output, but I don't bother to check.
    // This will give the most consistent and clearn output
    const std::string m_error_color_start = "";
    const std::string m_warning_color_start = "";
    const std::string m_color_end = "";
#else
    const std::string m_error_color_start = "\033[1;31m";
    const std::string m_warning_color_start = "\033[1;33m";
    const std::string m_color_end = "\033[0m";
#endif
};

#endif /* REPORTER_HPP_ */
