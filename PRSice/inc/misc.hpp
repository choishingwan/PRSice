//
//  misc.hpp
//  plink
//
//  Created by Shing Wan Choi on 18/08/2016.
//  Copyright © 2016 Shing Wan Choi. All rights reserved.
//

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

namespace misc{
    // my own codes
	inline bool to_bool(const std::string &str){
		if(str.compare("true")==0 || str.compare("T")==0 || str.compare("True")==0 ||
				str.compare("1")==0)
			return true;
		else if(str.compare("false")==0 || str.compare("F")==0 || str.compare("False")==0 ||
				str.compare("0")==0)
			return false;
		else{
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
    T convert(const std::string& str){
        std::istringstream iss(str);
        T obj;
        
        iss >> std::ws >> obj >> std::ws;
        
        if(!iss.eof())
            throw std::runtime_error("Unable to convert the input");
        
        return obj;
    }
    // trim from start (in place)
    inline void ltrim(std::string &s){
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    };
    // trim from end (in place)
    inline void rtrim(std::string &s){
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    };
    // trim from both ends (in place)
    inline void trim(std::string &s){
        ltrim(s);
        rtrim(s);
    };
    // trim from start (copying)
    inline std::string ltrimmed(std::string s){
        ltrim(s);
        return s;
    };
    // trim from end (copying)
    inline std::string rtrimmed(std::string s){
        rtrim(s);
        return s;
    };
    // trim from both ends (copying)
    inline std::string trimmed(std::string s){
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

}
#endif /* misc_hpp */
