/*
 * Reporter.hpp
 *
 *  Created on: 1 Nov 2017
 *      Author: shingwan
 */

#ifndef REPORTER_HPP_
#define REPORTER_HPP_

#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

class Reporter
{
public:
    Reporter(){

    };
    Reporter(const std::string& log_name, size_t width = 80) : m_width(width)
    {
        m_log_file.open(log_name.c_str());
        if (!m_log_file.is_open()) {
            std::string error_message =
                "Error: " + log_name + " cannot be open";
            throw std::runtime_error(error_message);
        }
    };
    void initiailize(const std::string& log_name, size_t width = 80)
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
    virtual ~Reporter()
    {
        if (m_log_file.is_open()) m_log_file.close();
    };
    void report(const std::string& input, bool wrap = true);

private:
    std::ofstream m_log_file;
    size_t m_width = 80;
};

#endif /* REPORTER_HPP_ */
