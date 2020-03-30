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


#include "misc.hpp"

namespace misc
{
std::istringstream Convertor::iss;

template <>
size_t Convertor::convert<size_t>(const std::string& str)
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
double Convertor::convert<double>(const std::string& str)
{
    iss.clear();
    iss.str(str);
    double obj;
    iss >> obj;
    if (!iss.eof() || iss.fail()
        || (std::fpclassify(obj) != FP_NORMAL
            && std::fpclassify(obj) != FP_ZERO)
        || errno == ERANGE)
    {
        std::cerr << "Checking: " << errno << "\t" << ERANGE << "\t" << str
                  << "\t" << std::fpclassify(obj) << "\t" << FP_NORMAL << "\t"
                  << FP_ZERO << std::endl;
        throw std::runtime_error("Unable to convert the input");
    }
    return obj;
}


double dnorm(double x, double mu, double sigma, bool log)
{
#ifdef IEEE_754
    if (isnan(x) != 0 || isnan(mu) != 0 || isnan(sigma) != 0)
        return x + mu + sigma;
#endif
    if (!(std::isfinite(sigma))) return 0.0;
    if (!(std::isfinite(x)) && mu == x) throw std::runtime_error("NA produced");
    if (sigma <= 0)
    {
        if (sigma < 0) throw std::runtime_error("Negative sigma not allowed");
        return (x == mu) ? std::numeric_limits<double>::infinity() : 0.0;
    }
    x = (x - mu) / sigma;
    if (!std::isfinite(x)) return 0.0;
    x = std::fabs(x);
    if (x >= 2 * std::sqrt(std::numeric_limits<double>::max())) return 0.0;
    if (log) return -(std::log(2 * M_PI) / 2 + 0.5 * x * x + std::log(sigma));
    if (x < 5)
        return (1 / std::sqrt(2 * M_PI)) * std::exp(-0.5 * x * x) / sigma;
    if (x > std::sqrt(-2 * std::log(2)
                      * (std::numeric_limits<double>::min_exponent + 1
                         - std::numeric_limits<double>::digits)))
        return 0.0;
    double x1 = std::ldexp(std::nearbyint(std::ldexp(x, 16)), -16);
    double x2 = x - x1;
    return (1 / std::sqrt(2 * M_PI)) / sigma
           * (exp(-0.5 * x1 * x1) * exp((-0.5 * x2 - x1) * x2));
}

double qnorm(double p, double mu, double sigma, bool lower_tail, bool log_p)
{
    double r, val;
#ifdef IEEE_754
    if (isnan(p) != 0 || isnan(mu) != 0 || isnan(sigma) != 0)
        return x + mu + sigma;
#endif
    if (log_p && p > 0.0)
        throw std::runtime_error("log p-value cannot be larger than 0");
    else if (!log_p && (p < 0.0 || p > 1.0))
        throw std::runtime_error("p-value out of bound");
    if (sigma < 0.0) throw std::runtime_error("Negative sigma not allowed");
    if (sigma == 0) return mu;
    double p_ = (log_p ? (lower_tail ? std::exp(p) : -expm1(p))
                       : (lower_tail ? (p) : (0.5 - (p) + 0.5)));
    double q = p_ - 0.5;
    if (fabs(q) <= .425) /* 0.075 <= p <= 0.925 */
    {
        r = .180625 - q * q;
        val = q
              * (((((((r * 2509.0809287301226727 + 33430.575583588128105) * r
                      + 67265.770927008700853)
                         * r
                     + 45921.953931549871457)
                        * r
                    + 13731.693765509461125)
                       * r
                   + 1971.5909503065514427)
                      * r
                  + 133.14166789178437745)
                     * r
                 + 3.387132872796366608)
              / (((((((r * 5226.495278852854561 + 28729.085735721942674) * r
                      + 39307.89580009271061)
                         * r
                     + 21213.794301586595867)
                        * r
                    + 5394.1960214247511077)
                       * r
                   + 687.1870074920579083)
                      * r
                  + 42.313330701600911252)
                     * r
                 + 1.);
    }
    else /* closer than 0.075 from {0,1} boundary */
    {
        /* r = min(p, 1-p) < 0.075 */
        if (q > 0)
            r = (log_p ? (lower_tail ? -expm1(p) : std::exp(p))
                       : (lower_tail ? (0.5 - (p) + 0.5) : (p))); /* 1-p */
        else
            r = p_; /* = R_DT_Iv(p) ^=  p */
        r = sqrt(-((log_p && ((lower_tail && q <= 0) || (!lower_tail && q > 0)))
                       ? p
                       : /* else */ log(r)));
        /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */
        if (r <= 5.) /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
        {
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 + .0227238449892691845833)
                            * r
                        + .24178072517745061177)
                           * r
                       + 1.27045825245236838258)
                          * r
                      + 3.64784832476320460504)
                         * r
                     + 5.7694972214606914055)
                        * r
                    + 4.6303378461565452959)
                       * r
                   + 1.42343711074968357734)
                  / (((((((r * 1.05075007164441684324e-9
                           + 5.475938084995344946e-4)
                              * r
                          + .0151986665636164571966)
                             * r
                         + .14810397642748007459)
                            * r
                        + .68976733498510000455)
                           * r
                       + 1.6763848301838038494)
                          * r
                      + 2.05319162663775882187)
                         * r
                     + 1.);
        }
        else /* very close to  0 or 1 */
        {
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7
                         + 2.71155556874348757815e-5)
                            * r
                        + .0012426609473880784386)
                           * r
                       + .026532189526576123093)
                          * r
                      + .29656057182850489123)
                         * r
                     + 1.7848265399172913358)
                        * r
                    + 5.4637849111641143699)
                       * r
                   + 6.6579046435011037772)
                  / (((((((r * 2.04426310338993978564e-15
                           + 1.4215117583164458887e-7)
                              * r
                          + 1.8463183175100546818e-5)
                             * r
                         + 7.868691311456132591e-4)
                            * r
                        + .0148753612908506148525)
                           * r
                       + .13692988092273580531)
                          * r
                      + .59983220655588793769)
                         * r
                     + 1.);
        }

        if (q < 0.0) val = -val;
        /* return (q >= 0.)? r : -r ;*/
    }
    return mu + sigma * val;
}
}
