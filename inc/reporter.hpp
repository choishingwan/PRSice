/*
 * Reporter.hpp
 *
 *  Created on: 1 Nov 2017
 *      Author: shingwan
 */

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
    Reporter(){
        m_error_prefix_size = m_error_prefix.size();
        m_warning_prefix_size = m_warning_prefix.size();
    }
    Reporter(const std::string& log_name, size_t width = 60) : m_width(width)
    {
        m_log_file.open(log_name.c_str());
        if (!m_log_file.is_open()) {
            std::string error_message =
                "Error: " + log_name + " cannot be open";
            throw std::runtime_error(error_message);
        }
        m_error_prefix_size = m_error_prefix.size();
        m_warning_prefix_size = m_warning_prefix.size();
    }
    void initiailize(const std::string& log_name, size_t width = 60)
    {
        m_width = width;
        if (m_log_file.is_open()) m_log_file.close();
        m_log_file.open(log_name.c_str());
        if (!m_log_file.is_open()) {
            std::string error_message =
                "Error: " + log_name + " cannot be open";
            throw std::runtime_error(error_message);
        }
    }
    virtual ~Reporter(){}
    void report(const std::string& input, bool wrap = true);

private:
    bool isNumeric(std::string s);
    std::ofstream m_log_file;
    std::string m_error_prefix = "Error:";
    std::string m_warning_prefix = "Warning:";
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
    std::string m_error_color_start = "";
    std::string m_warning_color_start = "";
    std::string m_color_end = "";
#else
    std::string m_error_color_start = "\033[1;31m";
    std::string m_warning_color_start = "\033[1;33m";
    std::string m_color_end = "\033[0m";
#endif
    size_t m_width = 60;
    size_t m_error_prefix_size;
    size_t m_warning_prefix_size;
};

#endif /* REPORTER_HPP_ */
