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
    std::vector<std::string> line;
    std::vector<std::string> check_list;
    split(paragraph, input.c_str(), '\n');
    std::string list_prefix;
    size_t space_pad = 0;
    for (auto&& message : paragraph)
    {
        check_list.clear();
        // go through each paragraph
        // red for error, yellow for warnings
        space_pad = 0;
        bool is_error =
            !message.compare(0, m_error_prefix_size, m_error_prefix);
        bool is_warning =
            !message.compare(0, m_warning_prefix_size, m_warning_prefix);
        bool is_list = false;
        if (!is_error && !is_warning)
        {
            // check if there are any list in this paragraph
            // list shouldn't be included in error or warning messages
            // list is assumed to be in paragraph format
            // (no in line list allowed)
            split(check_list, message.c_str(), ')');
            if (check_list.size() > 1)
            {
                if (misc::isNumeric(check_list.front()))
                {
                    is_list = true;
                    list_prefix = check_list.front() + ")";
                    space_pad = list_prefix.size() + 1;
                }
            }
        }
        if (is_error)
        {
            std::cerr << m_error_color_start;
            space_pad = m_error_prefix_size + 1;
        }
        else if (is_warning)
        {
            std::cerr << m_warning_color_start;
            space_pad = m_warning_prefix_size + 1;
        }
        if (!wrap)
        {
            std::cerr << message;
            m_log_file << message;
        }
        else
        {
            // now split the line into word and print it
            // out nicely
            split(line, message.c_str(), ' ');
            size_t cur_length = 0;
            for (auto&& word : line)
            {
                if (word.length() + cur_length >= m_width)
                {
                    // word is too long, so display on the next line
                    cur_length = word.length() + 1 + space_pad;
                    std::cerr << '\n';
                    m_log_file << '\n';
                    std::cerr << std::string(space_pad, ' ');
                    m_log_file << std::string(space_pad, ' ');
                }
                // line still have space. Add the word directly
                else
                    cur_length += word.length() + 1;
                std::cerr << word << " ";
                m_log_file << word << " ";
            }
        }
        if (is_error || is_warning) std::cerr << m_color_end;
        line.clear();
        std::cerr << '\n';
        m_log_file << '\n';
    }
    std::cerr << '\n' << std::flush;
    m_log_file << '\n' << std::flush;
}
