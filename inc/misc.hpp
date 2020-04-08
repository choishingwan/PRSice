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

#pragma once

#include <assert.h>
#include <stdexcept>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <algorithm>
#include <cmath>
#include <gzstream.h>
#include <iostream>
#include <limits>
#include <math.h>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>
#if defined __APPLE__
#include <mach/mach.h>
#include <mach/mach_host.h>
#include <mach/mach_init.h>
#include <mach/mach_types.h>
#include <mach/vm_statistics.h>
#include <sys/sysctl.h>
#elif defined _WIN32
#include <windows.h>
// psapi must go after windows, or will generate error
#include "psapi.h"
#else
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <sys/param.h>
#include <sys/sysctl.h>
#include <sys/types.h>
#include <unistd.h>
#endif
#if defined(_WIN32)
#include <windows.h>

#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) \
    || (defined(__APPLE__) && defined(__MACH__))
#include <sys/resource.h>
#include <unistd.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) \
    || (defined(__sun__) || defined(__sun)     \
        || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) \
    || defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif
#if defined(_WIN32)
#elif defined(__unix__) || defined(__unix) || defined(unix) \
    || (defined(__APPLE__) && defined(__MACH__))
#include <sys/param.h>
#include <sys/types.h>
#include <unistd.h>
#if defined(BSD)
#include <sys/sysctl.h>
#endif

#else
#error "Unable to define getMemorySize( ) for an unknown OS."
#endif
#define BIGSTACK_MIN_MB 64
#define BIGSTACK_DEFAULT_MB 2048


namespace misc
{
// my own codes

template <class T>
class vec2d
{
public:
    vec2d() {}
    vec2d(size_t row, size_t col, T def)
    {
        if (row == 0 || col == 0)
        { throw std::invalid_argument("Dimension of 2d vector must be >0"); }
        m_storage.resize(row * col, def);
        m_row = row;
        m_col = col;
    }
    vec2d(size_t row, size_t col)
    {
        if (row == 0 || col == 0)
        { throw std::invalid_argument("Dimension of 2d vector must be >0"); }
        m_storage.resize(row * col);
        m_row = row;
        m_col = col;
    }
    T operator()(size_t row, size_t col) const
    {
        if (row > m_row || col > m_col)
            throw std::out_of_range("2d vector out of range!");
        return m_storage[row * m_col + col];
    }
    T& operator()(size_t row, size_t col)
    {
        if (row > m_row || col > m_col)
            throw std::out_of_range("2d vector out of range!");
        return m_storage[row * m_col + col];
    }
    void clear() { m_storage.clear(); }
    size_t rows() const { return m_row; }
    size_t cols() const { return m_col; }

private:
    size_t m_row = 0;
    size_t m_col = 0;
    std::vector<T> m_storage;
};


inline bool to_bool(const std::string& input)
{
    std::string str = input;
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    if (str.compare("T") == 0 || str.compare("TRUE") == 0)
        return true;
    else if (str.compare("F") == 0 || str.compare("FALSE") == 0)
        return false;
    else
    {
        std::string error_message = "Input is not True/False values: " + str;
        throw std::runtime_error(error_message);
    }
}

template <typename T>
inline bool within_bound(const T& input, const T& low_bound, const T& up_bound)
{
    assert(low_bound <= up_bound);
    return !(input < low_bound || input > up_bound);
}

// TODO: Delete this, doesn't seems to give robust answer
inline size_t current_ram_usage() { return 0; }
// TODO: Delete this, doesn't seems to give robust answer
inline size_t total_ram_available()
{
#ifdef __APPLE__
    int32_t mib[2];
    size_t sztmp;
#endif
    unsigned char* bigstack_ua = nullptr; // ua = unaligned
    int64_t llxx;
    intptr_t default_alloc_mb;
    intptr_t malloc_size_mb = 0;
#ifdef __APPLE__
    mib[0] = CTL_HW;
    mib[1] = HW_MEMSIZE;
    llxx = 0;

    sztmp = sizeof(int64_t);
    sysctl(mib, 2, &llxx, &sztmp, nullptr, 0);
    llxx /= 1048576;
#else
#ifdef _WIN32
    MEMORYSTATUSEX memstatus;
    memstatus.dwLength = sizeof(memstatus);
    GlobalMemoryStatusEx(&memstatus);
    llxx = memstatus.ullTotalPhys / 1048576;
#else
    llxx = ((uint64_t) sysconf(_SC_PHYS_PAGES))
           * ((size_t) sysconf(_SC_PAGESIZE)) / 1048576;
#endif
#endif
    if (!llxx) { default_alloc_mb = BIGSTACK_DEFAULT_MB; }
    else if (llxx < (BIGSTACK_MIN_MB * 2))
    {
        default_alloc_mb = BIGSTACK_MIN_MB;
    }
    else
    {
        default_alloc_mb = llxx / 2;
    }
    if (!malloc_size_mb) { malloc_size_mb = default_alloc_mb; }
    else if (malloc_size_mb < BIGSTACK_MIN_MB)
    {
        malloc_size_mb = BIGSTACK_MIN_MB;
    }
    std::string message = "";
#ifndef __LP64__
    if (malloc_size_mb > 2047) { malloc_size_mb = 2047; }
#endif
    bigstack_ua =
        (unsigned char*) malloc(malloc_size_mb * 1048576 * sizeof(char));
    // if fail, return nullptr which will then get into the while loop
    while (!bigstack_ua)
    {
        malloc_size_mb = (malloc_size_mb * 3) / 4;
        if (malloc_size_mb < BIGSTACK_MIN_MB)
        { malloc_size_mb = BIGSTACK_MIN_MB; }
        bigstack_ua =
            (unsigned char*) malloc(malloc_size_mb * 1048576 * sizeof(char));
        if (bigstack_ua) {}
        else if (malloc_size_mb == BIGSTACK_MIN_MB)
        {
            throw std::runtime_error("Failed to allocate required memory");
        }
    }
    free(bigstack_ua);
    bigstack_ua = nullptr;
    return malloc_size_mb * 1024 * 1024;
}
// function from John D.Cook
// https://www.johndcook.com/blog/standard_deviation/
class RunningStat
{
public:
    RunningStat() {}
    void clear()
    {
        n = 0;
        M1 = M2 = M3 = M4 = 0.0;
    }
    void push(double x)
    {
        double delta, delta_n, delta_n2, term1;

        size_t n1 = n;
        n++;
        delta = x - M1;
        assert(n > 0);
        delta_n = delta / n;
        delta_n2 = delta_n * delta_n;
        term1 = delta * delta_n * n1;
        M1 += delta_n;
        M4 += term1 * delta_n2 * (n * n - 3 * n + 3) + 6 * delta_n2 * M2
              - 4 * delta_n * M3;
        M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
        M2 += term1;
    }
    size_t get_n() const { return n; }

    double mean() const { return M1; }

    double var() const { return M2 / (static_cast<double>(n) - 1.0); }

    double sd() const { return sqrt(var()); }

private:
    size_t n = 0;
    double M1 = 0, M2 = 0, M3 = 0, M4 = 0;
};


// Functions from R
double dnorm(double x, double mu = 0.0, double sigma = 1.0, bool log = false);
double qnorm(double p, double mu = 0.0, double sigma = 1.0,
             bool lower_tail = true, bool log_p = false);

// codes from stackoverflow
inline std::vector<std::string> split(const std::string& seq,
                                      const std::string& separators = "\t ")
{
    std::size_t prev = 0, pos;
    std::vector<std::string> result;
    while ((pos = seq.find_first_of(separators, prev)) != std::string::npos)
    {
        if (pos > prev) result.emplace_back(seq.substr(prev, pos - prev));
        prev = pos + 1;
    }
    if (prev < seq.length())
        result.emplace_back(seq.substr(prev, std::string::npos));
    return result;
}

inline std::vector<std::string_view> tokenize(std::string_view str,
                                              std::string delims = "\t ")
{
    std::vector<std::string_view> output;
    // output.reserve(str.size() / 2);
    for (auto first = str.data(), second = str.data(),
              last = first + str.size();
         second != last && first != last; first = second + 1)
    {
        second = std::find_first_of(first, last, std::cbegin(delims),
                                    std::cend(delims));
        if (first != second) output.emplace_back(first, second - first);
    }
    return output;
}

inline std::vector<std::string_view> split(std::string_view str,
                                           std::string delims = " ")
{
    std::vector<std::string_view> output;
    // output.reserve(str.size() / 2);
    for (auto first = str.data(), second = str.data(),
              last = first + str.size();
         second != last && first != last; first = second + 1)
    {
        second = std::find_first_of(first, last, std::cbegin(delims),
                                    std::cend(delims));
        if (first != second) output.emplace_back(first, second - first);
    }
    return output;
}

inline void split(std::vector<std::string>& result, const std::string& seq,
                  const std::string& separators = "\t ")
{
    std::size_t prev = 0, pos, idx = 0;
    const size_t init_size = result.size();
    // assuming we have the same size
    // result.clear();
    while ((pos = seq.find_first_of(separators, prev)) != std::string::npos)
    {
        if (pos > prev)
        {
            if (idx >= init_size)
            { result.emplace_back(seq.substr(prev, pos - prev)); }
            else
            {
                result[idx] = seq.substr(prev, pos - prev);
            }
            ++idx;
        }
        prev = pos + 1;
    }
    if (prev < seq.length())
    {
        if (idx > init_size)
        { result.emplace_back(seq.substr(prev, std::string::npos)); }
        else
        {
            result[idx] = seq.substr(prev, std::string::npos);
        }
        ++idx;
    }
    if (idx < init_size) { result.resize(idx); }
}

class Convertor
{
public:
    template <typename T>
    static T convert(const std::string& str, bool verbose = false)
    {
        iss.clear();
        iss.str(str);
        T obj;
        iss >> obj;
        if (!iss.eof() || iss.fail())
        { throw std::runtime_error("Unable to convert the input"); }
        return obj;
    }


private:
    static std::istringstream iss;
};

template <typename T>
inline T convert(const std::string& str, bool verbose = false)
{
    std::istringstream iss(str);
    T obj;
    iss >> obj;

    if (!iss.eof() || iss.fail())
    { throw std::runtime_error("Unable to convert the input"); }
    return obj;
}

template <>
inline size_t Convertor::convert<size_t>(const std::string& str, bool verbose)
{
    iss.clear();
    iss.str(str);
    size_t obj;
    iss >> obj;
    if (!iss.eof() || iss.fail())
    { throw std::runtime_error("Unable to convert the input"); }
    else if (static_cast<int>(obj) < 0)
    {
        throw std::runtime_error(
            "Error: Negative input for a positive "
            "variable, or you have a very large integer, e.g. larger than "
            + std::to_string(std::numeric_limits<int>::max()));
    }
    return obj;
}

template <>
inline double Convertor::convert<double>(const std::string& str, bool verbose)
{
    errno = 0;
    iss.clear();
    iss.str(str);
    double obj;
    iss >> obj;
    if (verbose)
    {
        std::cerr << "check: " << obj << "\t" << str << "\t"
                  << std::fpclassify(obj) << "\t" << FP_NORMAL << "\t"
                  << FP_ZERO << "\t" << errno << "\t" << ERANGE << std::endl;
    }
    if (!iss.eof() || iss.fail()
        || (std::fpclassify(obj) != FP_NORMAL
            && std::fpclassify(obj) != FP_ZERO)
        || errno == ERANGE)
    { throw std::runtime_error("Unable to convert the input"); }

    return obj;
}

template <typename T>
inline std::string to_string(T value)
{
    std::stringstream out;
    out << value;
    return out.str();
}

// NOTE: Didn't work for non-ASCII characters
inline void to_upper(std::string& str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}
inline void to_lower(std::string& str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
}

inline void to_upper(const std::string& input, std::string& out)
{
    out.resize(input.size());
    std::transform(input.begin(), input.end(), out.begin(), ::toupper);
}
inline void to_lower(const std::string& input, std::string& out)
{
    out.resize(input.size());
    std::transform(input.begin(), input.end(), out.begin(), ::tolower);
}
template <class T>
inline bool overflow(const T a, const T b)
{
    if (a == 0 || b == 0) return false;
    T result = a * b;
    return !(a == result / b);
}


// trim functions from https://stackoverflow.com/a/217605
// trim from start (in place)
inline std::string_view ltrimmed(std::string_view s)
{
    s.remove_prefix(std::distance(
        s.cbegin(), std::find_if(s.cbegin(), s.cend(),
                                 [](int ch) { return std::isgraph(ch); })));

    return s;
}

inline std::string_view rtrimmed(std::string_view s)
{
    s.remove_suffix(std::distance(
        s.crbegin(), std::find_if(s.crbegin(), s.crend(),
                                  [](int ch) { return std::isgraph(ch); })));
    return s;
}

inline std::string_view trimmed(std::string_view s)
{
    return ltrimmed(rtrimmed(s));
}
inline void ltrim(std::string_view& s)
{
    s.remove_prefix(std::distance(
        s.cbegin(), std::find_if(s.cbegin(), s.cend(),
                                 [](int ch) { return std::isgraph(ch); })));
}

inline void rtrim(std::string_view& s)
{
    s.remove_suffix(std::distance(
        s.crbegin(), std::find_if(s.crbegin(), s.crend(),
                                  [](int ch) { return std::isgraph(ch); })));
}

inline void trim(std::string_view& s)
{
    rtrim(s);
    ltrim(s);
}


inline void ltrim(std::string& s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                    [](int ch) { return std::isgraph(ch); }));
};

// trim from end (in place)
inline void rtrim(std::string& s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         [](int ch) { return std::isgraph(ch); })
                .base(),
            s.end());
};
// trim from both ends (in place)
inline void trim(std::string& s)
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
template <class T>
inline T base_name(T const& path, T const& delims = "/\\")
{
    return path.substr(path.find_last_of(delims) + 1);
}

template <class T>
inline T remove_extension(T const& filename)
{
    typename T::size_type const p(filename.find_last_of('.'));
    return p > 0 && p != T::npos ? filename.substr(0, p) : filename;
}

inline void replace_substring(std::string& s, const std::string& search,
                              const std::string& replace)
{
    if (search.empty()) return;
    for (size_t pos = 0;; pos += replace.length())
    {
        // Locate the substring to replace
        pos = s.find(search, pos);
        if (pos == std::string::npos) break;
        // Replace by erasing and inserting
        s.erase(pos, search.length());
        s.insert(pos, replace);
    }
}


// From the amazing Chris Chang (PLINK 2)

#ifndef DBL_MAX
#define DBL_MAX 1.7976931348623157e308
#endif

#define MAXV(aa, bb) (((bb) > (aa)) ? (bb) : (aa))
#define MINV(aa, bb) (((bb) < (aa)) ? (bb) : (aa))
const double kLentzFpmin = 1.0e-30;
const double kBigEpsilon = 0.000000476837158203125;
const double kSqrt2 = 1.4142135623730951;
const double kSqrtPi = 1.7724538509055159;
const double kLanczosG = 5.581;
const double log_min_value = -708.0;
const double log_max_value = 709.0;
const double kRecipE = 0.36787944117144233;
const double kLanczosSumExpgNumer[6] = {32.812445410297834, 32.123889414443320,
                                        12.580347294552161, 2.4634444783532414,
                                        0.2412010548258800, 0.0094469677045392};
const double kLanczosSumExpgDenom[6] = {0, 24, 50, 35, 10, 1};
const double kSmallEpsilon = 0.00000000000005684341886080801486968994140625;
const double kTemmeC0[7] = {-0.333333333,  0.0833333333,   -0.0148148148,
                            0.00115740741, 0.000352733686, -0.000178755144,
                            0.391926318e-4};
const double kTemmeC1[5] = {-0.00185185185, -0.00347222222, 0.00264550265,
                            -0.000990226337, 0.000205761317};
const double kTemmeC2[3] = {0.00413359788, -0.00268132716, 0.000771604938};
const double kFactorials[30] = {1.0,
                                1.0,
                                2.0,
                                6.0,
                                24.0,
                                120.0,
                                720.0,
                                5040.0,
                                40320.0,
                                362880.0,
                                3628800.0,
                                39916800.0,
                                479001600.0,
                                6227020800.0,
                                87178291200.0,
                                1307674368000.0,
                                20922789888000.0,
                                355687428096000.0,
                                6402373705728000.0,
                                121645100408832000.0,
                                0.243290200817664e19,
                                0.5109094217170944e20,
                                0.112400072777760768e22,
                                0.2585201673888497664e23,
                                0.62044840173323943936e24,
                                0.15511210043330985984e26,
                                0.403291461126605635584e27,
                                0.10888869450418352160768e29,
                                0.304888344611713860501504e30,
                                0.8841761993739701954543616e31};


inline uint32_t realnum(double dd)
{
    return (dd == dd) && (dd != INFINITY) && (dd != -INFINITY);
};

inline double betacf_slow(double aa, double bb, double xx)
{
    double qab = aa + bb;
    double qap = aa + 1.0;
    double qam = aa - 1.0;
    double cc = 1.0;
    double dd = 1.0 - qab * xx / qap;
    if (fabs(dd) < kLentzFpmin) { dd = kLentzFpmin; }
    dd = 1.0 / dd;
    double hh = dd;
    // evaluate 1 / (1 + d_1 / (1 + d_2 / (1 + d_3 / (...))))
    for (double mm = 1.0; mm <= 100.0; mm += 1.0)
    {
        double m2 = 2 * mm;

        // d_{2m}
        double tmp_aa = mm * (bb - mm) * xx / ((qam + m2) * (aa + m2));

        dd = 1.0 + tmp_aa * dd;
        if (fabs(dd) < kLentzFpmin) { dd = kLentzFpmin; }
        cc = 1.0 + tmp_aa / cc;
        if (fabs(cc) < kLentzFpmin) { cc = kLentzFpmin; }
        dd = 1.0 / dd;
        hh *= dd * cc;

        // d_{2m+1}
        tmp_aa = -(aa + mm) * (qab + mm) * xx / ((aa + m2) * (qap + m2));

        dd = 1.0 + tmp_aa * dd;
        if (fabs(dd) < kLentzFpmin) { dd = kLentzFpmin; }
        cc = 1.0 + tmp_aa / cc;
        if (fabs(cc) < kLentzFpmin) { cc = kLentzFpmin; }
        dd = 1.0 / dd;
        double del = dd * cc;
        hh *= del;
        if (fabs(del - 1.0) < 3.0e-7) { return hh; }
    }
    // don't detect failure for now
    return hh;
}

inline double betai_slow(double aa, double bb, double xx)
{
    if ((xx < 0.0) || (xx > 1.0)) { return -9; }
    uint32_t do_invert = (xx * (aa + bb + 2.0)) >= (aa + 1.0);
    if ((xx == 0.0) || (xx == 1.0)) { return (double) ((int32_t) do_invert); }
    // this is very expensive
    double bt = exp(lgamma(aa + bb) - lgamma(aa) - lgamma(bb) + aa * log(xx)
                    + bb * log(1.0 - xx));

    if (!do_invert) { return bt * betacf_slow(aa, bb, xx) / aa; }
    return 1.0 - bt * betacf_slow(bb, aa, 1.0 - xx) / bb;
}

inline double calc_tprob(double tt, double df)
{
    // must be thread-safe, so dcdflib won't cut it.
    // move this to plink2_stats once it's ready (and probably just eliminate
    // dcdflib entirely)
    if (!realnum(tt)) { return -9; }
    return betai_slow(df * 0.5, 0.5, df / (df + tt * tt));
}


inline double erfc_fast(double zz)
{
    const double tt = 1.0 / (1.0 + 0.5 * zz);
    const double tau =
        tt
        * exp(((((((((0.17087277 * tt - 0.82215223) * tt + 1.48851587) * tt
                    - 1.13520398)
                       * tt
                   + 0.27886807)
                      * tt
                  - 0.18628806)
                     * tt
                 + 0.09678418)
                    * tt
                + 0.37409196)
                   * tt
               + 1.00002368)
                  * tt
              - 1.26551223 - zz * zz);
    return tau;
}

inline double upper_gamma_fraction(double a1, double z1)
{
    // evaluate a_1 / (b_1 + (a_2 / (b_2 + (a_3 / (b_3 + ...)))))
    // see Boost continued_fraction_a(), upper_incomplete_gamma_fract
    double cur_b = z1 - a1 + 3;

    double hh = cur_b;
    const double a0 = a1 - 1.0;
    if (fabs(hh) < kLentzFpmin) { hh = kLentzFpmin; }
    double cc = hh;
    double dd = 0.0;
    for (double kk = 2.0; kk <= 100.0; kk += 1.0)
    {
        const double cur_a = kk * (a1 - kk);
        cur_b += 2.0;
        dd = cur_b + cur_a * dd;
        if (fabs(dd) < kLentzFpmin) { dd = kLentzFpmin; }
        cc = cur_b + cur_a / cc;
        if (fabs(cc) < kLentzFpmin) { cc = kLentzFpmin; }
        dd = 1.0 / dd;
        const double delta = cc * dd;
        hh *= delta;
        if (fabs(delta - 1.0) < 3.0e-7) { break; }
    }
    const double cont_frac = a0 / hh;
    return 1 / (z1 - a1 + 1 + cont_frac);
}

inline double small_gamma2_series(double aa, double xx, double init_value)
{
    double apn = aa + 1;
    const double negx = -xx;
    double nn = 1;
    double result = negx;
    double total = init_value;
    double rr;
    do
    {
        rr = result / apn;
        result *= negx;
        nn += 1.0;
        result /= nn;
        apn += 1;
        total += rr;
    } while (fabs(rr) > (kBigEpsilon * kBigEpsilon));
    return total;
}

inline double tgamma_small_upper_part_df1(double xx, uint32_t invert,
                                          double* p_derivative, double* pgam)
{
    // x < 1.1
    // df == 1, a == 0.5
    double result = 0.5 * kSqrtPi - 1.0;
    *pgam = (result + 1) * 2;
    double pp = sqrt(xx) - 1.0; // no point in using powm1() with ^0.5
    result -= pp;
    result *= 2;
    pp += 1;
    if (p_derivative) { *p_derivative = pp / ((*pgam) * exp(xx)); }
    const double init_value = invert ? (*pgam) : 0;
    result = -pp * small_gamma2_series(0.5, xx, (init_value - result) / pp);
    if (invert) { result = -result; }
    return result;
}

inline double lower_gamma_series(double aa, double zz, double init_value)
{
    // z must not be much larger than a
    double result = 1;
    double total = init_value;
    double rr;
    do
    {
        rr = result;
        aa += 1.0;
        result *= zz / aa;
        total += rr;
    } while (std::abs(rr) > (kBigEpsilon * kBigEpsilon));
    return total;
}

inline double finite_half_gamma_q(double aa, double xx, double* p_derivative)
{
    // a is in {0.5, 1.5, ..., 29.5}; max(0.2, a-1) < x < log_max_value
    const double sqrt_x = sqrt(xx);
    double ee = erfc_fast(sqrt_x);
    if ((ee != 0) && (aa > 1))
    {
        double term = exp(-xx) / (kSqrtPi * sqrt_x);
        term *= xx * 2;
        double sum = term;
        for (double nn = 1.5; nn < aa; nn += 1.0)
        {
            term /= nn;
            term *= xx;
            sum += term;
        }
        ee += sum;
        if (p_derivative) { *p_derivative = 0; }
    }
    else if (p_derivative)
    {
        *p_derivative = sqrt_x * exp(-xx) * (1.0 / kSqrtPi);
    }
    return ee;
}

inline double finite_gamma_q(uint32_t aa, double xx, double* p_derivative)
{
    // a is a positive integer < 30; max(0.6, a-1) < x < log_max_value
    // (e^{-x})(1 + x + x^2/2 + x^3/3! + x^4/4! + ... + x^{a-1}/(a-1)!)
    const double ee = exp(-xx);
    if (ee == 0.0) { return 0; }
    double sum = ee;
    double term = sum;
    for (uint32_t nn = 1; nn < aa; ++nn)
    {
        term /= (double) ((int32_t) nn);
        term *= xx;
        sum += term;
    }
    if (p_derivative)
    { *p_derivative = ee * pow(xx, (int32_t) aa) / kFactorials[aa - 1]; }
    return sum;
}

inline double lanczos_sum_expg_scaled_recip(double zz)
{
    double s1;
    double s2;
    if (zz <= 1)
    {
        s1 = kLanczosSumExpgNumer[5];
        s2 = kLanczosSumExpgDenom[5];
        for (int32_t ii = 4; ii >= 0; --ii)
        {
            s1 *= zz;
            s2 *= zz;
            s1 += kLanczosSumExpgNumer[(uint32_t) ii];
            s2 += kLanczosSumExpgDenom[(uint32_t) ii];
        }
    }
    else
    {
        zz = 1 / zz;
        s1 = kLanczosSumExpgNumer[0];
        s2 = kLanczosSumExpgDenom[0];
        for (uint32_t uii = 1; uii < 6; ++uii)
        {
            s1 *= zz;
            s2 *= zz;
            s1 += kLanczosSumExpgNumer[uii];
            s2 += kLanczosSumExpgDenom[uii];
        }
    }
    // may as well flip this
    return s2 / s1;
}

inline double log1pmx(double xx)
{
    // log(1+x) - x
    // assumes abs(xx) < 0.95
    const double aa = fabs(xx);
    if (aa < (kBigEpsilon / kSqrt2))
    { // 2^{-21.5}
        return -xx * xx * 0.5;
    }
    double kk = 1.0; // skip first term of usual log(1+x) series
    const double m_mult = -xx;
    double m_prod = xx;
    double total = 0.0;
    double rr;
    do
    {
        m_prod *= m_mult;
        kk += 1.0;
        rr = m_prod / kk;
        total += rr;
        // todo: tune these epsilons, but let's wait until we know all of the
        // callers of these functions
    } while (fabs(rr) > (kBigEpsilon * kBigEpsilon));
    return total;
}

inline double regularized_gamma_prefix(double aa, double zz)
{
    // assumes a == 0.5 if a < 1.  assumes z > 0.
    // we are fine with float-level precision, so lanczos_n=6, kLanczosG=5.581
    if (aa < 1) { return sqrt(zz) * exp(-zz) * (1.0 / kSqrtPi); }
    const double agh = aa + kLanczosG - 0.5;
    const double agh_recip = 1.0 / agh;
    const double dd = ((zz - aa) - (kLanczosG - 0.5)) * agh_recip;
    double prefix;
    if ((fabs(dd * dd * aa) <= 100) && (aa > 150))
    {
        // abs(dd) < sqrt(2/3) < 0.95
        prefix = aa * log1pmx(dd) + zz * (0.5 - kLanczosG) * agh_recip;
        prefix = exp(prefix);
    }
    else
    {
        const double alz = aa * log(zz * agh_recip);
        const double amz = aa - zz;
        const double cur_minv = MINV(alz, amz);
        if ((cur_minv <= log_min_value) || (MAXV(alz, amz) >= log_max_value))
        {
            const double amza = amz / aa;
            double sq;
            if ((cur_minv > 2 * log_min_value)
                && (MAXV(alz, amz) < 2 * log_max_value))
            {
                sq = pow(zz * agh_recip, aa * 0.5) * exp(amz * 0.5);
                prefix = sq * sq;
            }
            else if ((cur_minv > 4 * log_min_value)
                     && (MAXV(alz, amz) < 4 * log_max_value) && (zz > aa))
            {
                sq = pow(zz * agh_recip, aa * 0.25) * exp(amz * 0.25);
                prefix = sq * sq;
                prefix *= prefix;
            }
            else if ((amza > log_min_value) && (amza < log_max_value))
            {
                prefix = pow((zz * exp(amza)) * agh_recip, aa);
            }
            else
            {
                prefix = exp(alz + amz);
            }
        }
        else
        {
            prefix = pow(zz * agh_recip, aa) * exp(amz);
        }
    }
    prefix *= sqrt(agh * kRecipE) * lanczos_sum_expg_scaled_recip(aa);
    return prefix;
}

inline double igamma_temme_large(double aa, double xx)
{
    // 24-bit precision is fine
    const double sigma = (xx - aa) / aa;
    // abs(sigma) < 0.4
    const double phi = -log1pmx(sigma);
    const double sqrt_a = sqrt(aa);
    const double sqrt_phi = sqrt(phi);
    const double yy = aa * phi;
    double zz = kSqrt2 * sqrt_phi;
    if (xx < aa) { zz = -zz; }
    double workspace[3];
    workspace[0] = (((((kTemmeC0[6] * zz + kTemmeC0[5]) * zz + kTemmeC0[4]) * zz
                      + kTemmeC0[3])
                         * zz
                     + kTemmeC0[2])
                        * zz
                    + kTemmeC0[1])
                       * zz
                   + kTemmeC0[0];
    workspace[1] = (((kTemmeC1[4] * zz + kTemmeC1[3]) * zz + kTemmeC1[2]) * zz
                    + kTemmeC1[1])
                       * zz
                   + kTemmeC1[0];
    workspace[2] = (kTemmeC2[2] * zz + kTemmeC2[1]) * zz + kTemmeC2[0];
    const double a_recip = 1 / aa;
    double result =
        (workspace[2] * a_recip + workspace[1]) * a_recip + workspace[0];
    result *= exp(-yy) / ((kSqrt2 * kSqrtPi) * sqrt_a);
    if (xx < aa) { result = -result; }
    result += erfc_fast(sqrt_a * sqrt_phi) * 0.5;
    return result;
}

inline double gamma_incomplete_imp2(uint32_t df, double xx, uint32_t invert,
                                    double* p_derivative)
{
    assert(df);
    assert(xx >= 0.0);
    const double aa = ((double) ((int32_t) df)) * 0.5;
    const uint32_t is_small_a =
        (df < 60) && (aa <= xx + 1) && (xx < log_max_value);
    uint32_t is_int = 0;
    uint32_t is_half_int = 0;
    if (is_small_a)
    {
        is_half_int = df % 2;
        is_int = !is_half_int;
    }
    uint32_t eval_method;
    if (is_int && (xx > 0.6))
    {
        invert = !invert;
        eval_method = 0;
    }
    else if (is_half_int && (xx > 0.2))
    {
        invert = !invert;
        eval_method = 1;
    }
    else if (xx < kSmallEpsilon)
    {
        // avoid computing log(0)
        // don't need more precision here, 6 digits is enough
        assert(!p_derivative);
        return 1.0;
    }
    else if (xx < 0.5)
    {
        // log(x) is negative
        // -0.4 / log(x) >= 0.5 (this is impossible for larger a)
        // -> -0.4 <= 0.5 * log(x)
        // -> -0.8 <= log(x)
        // -> e^{-0.8} <= x
        eval_method = 2 + ((df == 1) && (xx >= 0.44932896411722156));
    }
    else if (xx < 1.1)
    {
        // x * 0.75 >= 0.5
        // x >= 2/3
        eval_method = 2 + ((df == 1) && (xx >= (2.0 / 3.0)));
    }
    else
    {
        const double x_minus_a = xx - aa;
        uint32_t use_temme = 0;
        if (aa > 20)
        {
            // sigma = abs((x - a) / a);
            // igamma_temme_large() assumes abs(sigma) < 0.95
            if (aa > 200)
            {
                // abs(sigma) < sqrt(20 / a) < 0.316...
                use_temme = (20 * aa > x_minus_a * x_minus_a);
            }
            else
            {
                // abs(sigma) < 0.4
                const double sigma_times_a = fabs(x_minus_a);
                use_temme = (sigma_times_a < 0.4 * aa);
            }
        }
        if (use_temme) { eval_method = 5; }
        else
        {
            // x - (1 / (3 * x)) < a
            // x * x - (1/3) < a * x
            // x * x - a * x < 1/3
            // x * (x - a) < 1/3
            if (xx * x_minus_a < (1.0 / 3.0)) { eval_method = 2; }
            else
            {
                eval_method = 4;
                invert = !invert;
            }
        }
    }
    double result;
    switch (eval_method)
    {
    case 0: result = finite_gamma_q(df / 2, xx, p_derivative); break;
    case 1:
        // previously used erfc, but that was up to ~3x as slow as dcdflib (e.g.
        // chiprob_p(2.706, 1) case).
        result = finite_half_gamma_q(aa, xx, p_derivative);
        if (p_derivative && (*p_derivative == 0))
        { *p_derivative = regularized_gamma_prefix(aa, xx); }
        break;
    case 2:
        result = regularized_gamma_prefix(aa, xx);
        if (p_derivative) { *p_derivative = result; }
        if (result != 0)
        {
            // uint32_t optimized_invert = 0;
            double init_value = 0;
            if (invert)
            {
                init_value = -aa / result;
                // optimized_invert = 1;
            }
            result *= lower_gamma_series(aa, xx, init_value) / aa;
            // if (optimized_invert) {
            if (invert)
            {
                invert = 0;
                result = -result;
            }
        }
        break;
    case 3:
    {
        invert = !invert;
        double gg;
        result = tgamma_small_upper_part_df1(xx, invert, p_derivative, &gg);
        invert = 0;
        result /= gg;
    }
    break;
    case 4:
        result = regularized_gamma_prefix(aa, xx);
        if (p_derivative) { *p_derivative = result; }
        if (result != 0) { result *= upper_gamma_fraction(aa, xx); }
        break;
    case 5:
        result = igamma_temme_large(aa, xx);
        if (xx >= aa) { invert = !invert; }
        if (p_derivative) { *p_derivative = regularized_gamma_prefix(aa, xx); }
    }
    if (result > 1) { result = 1; }
    if (invert) { result = 1 - result; }
    if (p_derivative)
    {
        if ((xx < 1) && (DBL_MAX * xx < (*p_derivative)))
        {
            *p_derivative = DBL_MAX / 2; // overflow; do we really need this?
        }
        else
        {
            *p_derivative /= xx;
        }
    }
    return result;
}

// inline double chiprob_p(double chisq, uint32_t df)
//{
//    // todo: figure out when we were depending on this to return -9, and
//    decide
//    // how to handle those situations now
//    return gamma_incomplete_imp2(df, chisq * 0.5, 1, nullptr);
//}

/*!
 * \brief Function to check if two double are equal from
 *        https://stackoverflow.com/a/4010279/1441789
 * \param a the first double
 * \param b the second double
 * \param error_factor level of error, should be of no concern to us at the
 *        moment
 * \return True if two double are equal
 */
inline bool logically_equal(double a, double b, double error_factor = 1.0)
{
    return ((a == b)
            || (std::abs(a - b) < std::abs(std::min(a, b))
                                      * std::numeric_limits<double>::epsilon()
                                      * error_factor));
}

inline bool is_gz_file(const std::string& name)
{
    const unsigned char gz_magic[2] = {0x1f, 0x8b};
    FILE* fp;
    if ((fp = fopen(name.c_str(), "rb")) == nullptr)
    { throw std::runtime_error("Error: Cannot open file - " + name); }
    unsigned char buf[2];
    if (fread(buf, 1, 2, fp) == 2)
    {
        if (buf[0] == gz_magic[0] && buf[1] == gz_magic[1]) { return true; }
        return false;
    }
    else
    {
        // can open the file, but can't read the magic number.
        return false;
    }
}

inline std::unique_ptr<std::istream> load_stream(const std::string& filepath,
                                                 bool& gz_input)
{
    gz_input = false;
    try
    {
        gz_input = misc::is_gz_file(filepath);
    }
    catch (const std::runtime_error& e)
    {
        throw std::runtime_error(e.what());
    }
    if (gz_input)
    {
        auto gz =
            std::make_unique<GZSTREAM_NAMESPACE::igzstream>(filepath.c_str());
        if (!gz->good())
        {
            throw std::runtime_error("Error: Cannot open file: " + filepath
                                     + " (gz) to read!\n");
        }
        return std::unique_ptr<std::istream>(*gz ? std::move(gz) : nullptr);
    }
    else
    {
        auto file = std::make_unique<std::ifstream>(filepath.c_str());
        if (!file->is_open())
        { throw std::runtime_error("Error: Cannot open file: " + filepath); }
        return std::unique_ptr<std::istream>(*file ? std::move(file) : nullptr);
    }
}

inline bool isNumeric(const std::string& s)
{
    try
    {
        convert<double>(s);
    }
    catch (...)
    {
        return false;
    }
    return true;
}

inline int string_to_int(const char* p)
{
    int x = 0;
    bool neg = false;
    if (*p == '-')
    {
        neg = true;
        ++p;
    }
    else if (*p == '+')
    {
        ++p;
    }
    else if (*p < '0' || *p > '9')
    {
        throw std::runtime_error("Error: Not an integer\n");
    }
    while (*p >= '0' && *p <= '9')
    {
        x = (x * 10) + (*p - '0');
        ++p;
    }
    if (neg) { x = -x; }
    return x;
}
inline size_t string_to_size_t(const char* p)
{
    size_t x = 0;
    if (*p == '-')
    {
        throw std::runtime_error(
            "Error: Negative value, cannot be assigned to unsigned integer\n");
    }
    else if (*p == '+')
    {
        ++p;
    }
    else if (*p < '0' || *p > '9')
    {
        throw std::runtime_error("Error: Not an integer\n");
    }
    while (*p >= '0' && *p <= '9')
    {
        x = (x * 10) + static_cast<size_t>(*p - '0');
        ++p;
    }
    return x;
}
// from https://stackoverflow.com/a/874160
inline bool hasEnding(const std::string& fullString, const std::string& ending)
{
    if (fullString.empty())
        throw std::runtime_error(
            "Error: Cannot look for ending of an empty string");
    else if (ending.empty())
        throw std::runtime_error(
            "Error: Undefined behaviour. Cannot look for empty ending in "
            "string");
    else if (fullString.length() >= ending.length())
    {
        return (fullString.compare(fullString.length() - ending.length(),
                                   ending.length(), ending)
                == 0);
    }
    else
    {
        return false;
    }
}

/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 * From https://stackoverflow.com/a/14927379
 */
inline size_t getPeakRSS()
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (size_t) info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) \
    || (defined(__sun__) || defined(__sun)     \
        || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ((fd = open("/proc/self/psinfo", O_RDONLY)) == -1)
        return (size_t) 0L; /* Can't open? */
    if (read(fd, &psinfo, sizeof(psinfo)) != sizeof(psinfo))
    {
        close(fd);
        return (size_t) 0L; /* Can't read? */
    }
    close(fd);
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) \
    || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage(RUSAGE_SELF, &rusage);
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t) rusage.ru_maxrss;
#else
    return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t) 0L; /* Unsupported. */
#endif
}


/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
inline size_t getCurrentRSS()
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo(GetCurrentProcess(), &info, sizeof(info));
    return (size_t) info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t) &info,
                  &infoCount)
        != KERN_SUCCESS)
        return (size_t) 0L; /* Can't access? */
    return (size_t) info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) \
    || defined(__gnu_linux__)
    /* Linux ---------------------------------------------------- */
    long rss = 0L;
    FILE* fp = NULL;
    if ((fp = fopen("/proc/self/statm", "r")) == NULL)
        return (size_t) 0L; /* Can't open? */
    if (fscanf(fp, "%*s%ld", &rss) != 1)
    {
        fclose(fp);
        return (size_t) 0L; /* Can't read? */
    }
    fclose(fp);
    return (size_t) rss * (size_t) sysconf(_SC_PAGESIZE);

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t) 0L; /* Unsupported. */
#endif
}


/**
 * Returns the size of physical memory (RAM) in bytes.
 */
inline size_t getMemorySize()
{
#if defined(_WIN32) && (defined(__CYGWIN__) || defined(__CYGWIN32__))
    /* Cygwin under Windows. ------------------------------------ */
    /* New 64-bit MEMORYSTATUSEX isn't available.  Use old 32.bit */
    MEMORYSTATUS status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatus(&status);
    return (size_t) status.dwTotalPhys;

#elif defined(_WIN32)
    /* Windows. ------------------------------------------------- */
    /* Use new 64-bit MEMORYSTATUSEX, not old 32-bit MEMORYSTATUS */
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    GlobalMemoryStatusEx(&status);
    return (size_t) status.ullTotalPhys;

#elif defined(__unix__) || defined(__unix) || defined(unix) \
    || (defined(__APPLE__) && defined(__MACH__))
    /* UNIX variants. ------------------------------------------- */
    /* Prefer sysctl() over sysconf() except sysctl() HW_REALMEM and HW_PHYSMEM
     */

#if defined(CTL_HW) && (defined(HW_MEMSIZE) || defined(HW_PHYSMEM64))
    int mib[2];
    mib[0] = CTL_HW;
#if defined(HW_MEMSIZE)
    mib[1] = HW_MEMSIZE; /* OSX. --------------------- */
#elif defined(HW_PHYSMEM64)
    mib[1] = HW_PHYSMEM64; /* NetBSD, OpenBSD. --------- */
#endif
    int64_t size = 0;    /* 64-bit */
    size_t len = sizeof(size);
    if (sysctl(mib, 2, &size, &len, NULL, 0) == 0) return (size_t) size;
    return 0L; /* Failed? */

#elif defined(_SC_AIX_REALMEM)
    /* AIX. ----------------------------------------------------- */
    return (size_t) sysconf(_SC_AIX_REALMEM) * (size_t) 1024L;

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)
    /* FreeBSD, Linux, OpenBSD, and Solaris. -------------------- */
    return (size_t) sysconf(_SC_PHYS_PAGES) * (size_t) sysconf(_SC_PAGESIZE);

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGE_SIZE)
    /* Legacy. -------------------------------------------------- */
    return (size_t) sysconf(_SC_PHYS_PAGES) * (size_t) sysconf(_SC_PAGE_SIZE);

#elif defined(CTL_HW) && (defined(HW_PHYSMEM) || defined(HW_REALMEM))
    /* DragonFly BSD, FreeBSD, NetBSD, OpenBSD, and OSX. -------- */
    int mib[2];
    mib[0] = CTL_HW;
#if defined(HW_REALMEM)
    mib[1] = HW_REALMEM;   /* FreeBSD. ----------------- */
#elif defined(HW_PYSMEM)
    mib[1] = HW_PHYSMEM; /* Others. ------------------ */
#endif
    unsigned int size = 0; /* 32-bit */
    size_t len = sizeof(size);
    if (sysctl(mib, 2, &size, &len, NULL, 0) == 0) return (size_t) size;
    return 0L; /* Failed? */
#endif /* sysctl and sysconf variants */

#else
    return 0L;          /* Unknown OS. */
#endif
}

inline unsigned long long remain_memory(const double& adjFactor = 0.8)
{
    return (misc::getMemorySize() * adjFactor - getCurrentRSS());
}
}
