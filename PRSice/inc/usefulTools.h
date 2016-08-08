// This file is part of SHREK, Snp HeRitability Estimate Kit
//
// Copyright (C) 2014-2015 Sam S.W. Choi <choishingwan@gmail.com>
//
// This Source Code Form is subject to the terms of the GNU General
// Public License v. 2.0. If a copy of the GPL was not distributed
// with this file, You can obtain one at
// https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html.

#ifndef USEFULTOOLS_H
#define	USEFULTOOLS_H

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <math.h>

/** \class usefulTools
 *  \brief a collection of small utility functions gathered in different places
 */
class usefulTools{
    public:
    /** Tokenize the string using based on the separator (by Thomas, my FYP supervisor) */
    static void tokenizer(const std::string seq, const std::string separators, std::vector<std::string>* result);
    static std::vector<std::string> split(const std::string seq, const std::string separators){
    /** dnorm from R */
	static double dnorm(const double x);
	/** qnorm from R */
	static double qnorm(const double p);
    /* Black magic got from http://stackoverflow.com/a/109025 */

    /** Functions for getting the sign of a value input, obtained from  <a href="http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c">stackoverflow</a>  */
	template <typename T> inline constexpr
	static int signum(T x, std::false_type is_signed) { return T(0) < x; }
    /** Functions for getting the sign of a value input, obtained from  <a href="http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c">stackoverflow</a>  */
	template <typename T> inline constexpr
	static int signum(T x, std::true_type is_signed) { return (T(0) < x) - (x < T(0)); }
    /** Functions for getting the sign of a value input, obtained from  <a href="http://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c">stackoverflow</a>  */
	template <typename T> inline constexpr
	static int signum(T x) { return signum(x, std::is_signed<T>()); }


    static void ltrim(std::string &s);
    // trim from end (in place)
    static void rtrim(std::string &s);
    // trim from both ends (in place)
    static void trim(std::string &s);
    // trim from start (copying)
    static std::string ltrimmed(std::string s);
    // trim from end (copying)
    static std::string rtrimmed(std::string s);
    // trim from both ends (copying)
    static std::string trimmed(std::string s);
};

#endif	/* USEFULTOOLS_H */

