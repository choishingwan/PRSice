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

#include "regression.hpp"
namespace Regression
{

// on purposely perform the copying of x
void linear_regression(const Eigen::VectorXd& y, const Eigen::MatrixXd& A,
                       double& p_value, double& r2, double& r2_adjust,
                       double& coeff, double& standard_error, size_t thread,
                       bool intercept)
{
    Eigen::setNbThreads(static_cast<int>(thread));
    // in more general cases, the following is needed (adding the intercept)
    // but here in PRSice, we always include the intercept, so we will skip
    // this to speed things up
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> z(A);
    Eigen::VectorXd beta = z.solve(y);
    Eigen::VectorXd fitted = A * beta;
    double rss = (A * beta - y).squaredNorm();
    double mss = (fitted.array() - fitted.mean()).pow(2).sum();
    r2 = mss / (mss + rss);
    Eigen::Index n = A.rows();
    Eigen::Index rank = z.rank();
    double rdf = static_cast<double>(n - rank);
    long df_int = intercept; // 0 false 1 true
    r2_adjust = 1.0 - (1.0 - r2) * (static_cast<double>(n - df_int) / rdf);
    double resvar = rss / rdf;
    Eigen::Index se_index = intercept;
    for (Eigen::Index ind = 0; ind < beta.rows(); ++ind)
    {
        if (z.colsPermutation().indices()(ind) == intercept)
        {
            se_index = ind;
            break;
        }
    }

    Eigen::MatrixXd R =
        z.matrixR().topLeftCorner(rank, rank).triangularView<Eigen::Upper>();
    Eigen::VectorXd se =
        ((R.transpose() * R).inverse().diagonal() * resvar).array().sqrt();
    // Remember, only the coefficient's order is wrong e.g. intercept at the end
    double tval = beta(intercept)
                  / se(se_index); // only interested in the one coefficient
    coeff = beta(intercept);
    standard_error = se(se_index);
    p_value = misc::calc_tprob(tval, n);
}

// This is an unsafe version of R's glm.fit
// unsafe as in I have skipped some of the checking
void glm(const Eigen::VectorXd& y, const Eigen::MatrixXd& x, double& p_value,
         double& r2, double& coeff, double& standard_error, size_t thread)
{
    Binomial family = Binomial();
    Eigen::setNbThreads(static_cast<int>(thread));
    GLM<Binomial> run_glm(x, y);
    run_glm.init_parms(family);
    run_glm.solve(family);
    r2 = run_glm.get_r2(family);
    run_glm.get_stat(1, p_value, coeff, standard_error);
}
}
