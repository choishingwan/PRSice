// This file is part of PRSice2.0, copyright (C) 2016-2017
// Shing Wan Choi, Jack Euesden, Cathryn M. Lewis, Paul F. O’Reilly
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


#ifndef misc_hpp
#define misc_hpp

#include <stdio.h>
#include <stdexcept>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
#include <limits>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iostream>

namespace misc
{
// my own codes

    template <class T> class vec2d
    {
        public:
            vec2d()
            {

            }
            vec2d(int row, int col, T def)
            {
                if(row==0 || col==0)
                {
                    throw std::invalid_argument("Dimension of 2d vector must be >0");
                }
                m_storage.resize(row*col, def);
                m_row = row;
                m_col = col;
            };
            vec2d(int row, int col)
            {
                if(row==0 || col==0)
                {
                    throw std::invalid_argument("Dimension of 2d vector must be >0");
                }
                m_storage.resize(row*col);
                m_row = row;
                m_col = col;
            };
            T operator()(int row, int col) const
            {
                if(row > m_row || col>m_col)
                {
                    throw std::out_of_range("2d vector out of range!");
                }
                return m_storage[row*m_col+col];
            }
            T &operator()(int row, int col)
            {
                if(row > m_row || col>m_col)
                {
                    throw std::out_of_range("2d vector out of range!");
                }
                return m_storage[row*m_col+col];
            }
            void clear()
            {
                m_storage.clear();
            }
            size_t rows() const { return m_row; };
            size_t cols() const { return m_col; };
        private:
            size_t m_row = 0;
            size_t m_col = 0;
            std::vector<T> m_storage;
    };

    inline bool to_bool(const std::string &str)
    {
        if(str.compare("true")==0 || str.compare("T")==0 || str.compare("True")==0 ||
                str.compare("1")==0)
            return true;
        else if(str.compare("false")==0 || str.compare("F")==0 || str.compare("False")==0 ||
                str.compare("0")==0)
            return false;
        else
        {
            std::string error_message="Input is not True/False values: "+str;
            throw std::runtime_error(error_message);
        }
    }


// Functions from R
    double dnorm(double x, double mu=0.0, double sigma=1.0, bool log=false);
    double qnorm(double p, double mu=0.0, double sigma=1.0, bool lower_tail=true, bool log_p=false);


// codes from stackoverflow
    std::vector<std::string> split(const std::string seq, const std::string separators="\t ");
    template <typename T> inline
    T convert(const std::string& str)
    {
        std::istringstream iss(str);
        T obj;

        iss >> std::ws >> obj >> std::ws;

        if(!iss.eof())
            throw std::runtime_error("Unable to convert the input");

        return obj;
    }
// trim from start (in place)
    inline void ltrim(std::string &s)
    {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    };
// trim from end (in place)
    inline void rtrim(std::string &s)
    {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    };
// trim from both ends (in place)
    inline void trim(std::string &s)
    {
        ltrim(s);
        rtrim(s);
    };
// trim from start (copying)
    inline std::string ltrimmed(std::string s)
    {
        ltrim(s);
        return s;
    };
// trim from end (copying)
    inline std::string rtrimmed(std::string s)
    {
        rtrim(s);
        return s;
    };
// trim from both ends (copying)
    inline std::string trimmed(std::string s)
    {
        trim(s);
        return s;
    };
// From http://stackoverflow.com/a/24386991/1441789
    template<class T>
    inline T base_name(T const & path, T const & delims = "/\\")
    {
        return path.substr(path.find_last_of(delims) + 1);
    }
    template<class T>
    inline T remove_extension(T const & filename)
    {
        typename T::size_type const p(filename.find_last_of('.'));
        return p > 0 && p != T::npos ? filename.substr(0, p) : filename;
    }

    inline void replace_substring( std::string &s, const std::string &search, const std::string &replace )
    {
        for( size_t pos = 0; ; pos += replace.length() )
        {
            // Locate the substring to replace
            pos = s.find( search, pos );
            if( pos == std::string::npos ) break;
            // Replace by erasing and inserting
            s.erase( pos, search.length() );
            s.insert( pos, replace );
        }
    }

}
#endif /* misc_hpp */
