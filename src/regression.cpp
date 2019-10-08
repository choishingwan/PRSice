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

// This is an unsafe version of R's glm.fit
// unsafe as in I have skipped some of the checking
void glm(const Eigen::VectorXd& y, const Eigen::MatrixXd& x, double& p_value,
         double& r2, double& coeff, double& standard_error, int thread)
{
    Binomial family = Binomial();
    Eigen::setNbThreads(thread);
    GLM<Binomial> run_glm(x, y, family);
    run_glm.init_parms();
    run_glm.solve();
    r2 = run_glm.get_r2();
    run_glm.get_stat(1, p_value, coeff, standard_error);
}

void fastLm(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, double& p_value,
            double& r2, double& r2_adjust, double& coeff,
            double& standard_error, int thread, bool intercept, int type)
{
    Eigen::setNbThreads(thread);
    Eigen::Index n = X.rows();
    if (n != y.rows()) { throw std::runtime_error("Error: Size mismatch"); }
    lm ans;
    switch (type)
    {
    case 0: ans = ColPivQR(X, y); break;
    case 1: ans = QR(X, y); break;
    case 2: ans = Llt(X, y); break;
    case 3: ans = Ldlt(X, y); break;
    case 4: ans = SVD(X, y); break;
    case 5: ans = SymmEigen(X, y); break;
    default: throw std::runtime_error("Error: Invalid regression type");
    }
    coeff = ans.coef()(1);
    Eigen::Index rank = ans.rank();
    Eigen::VectorXd resid = y - ans.fitted();
    Eigen::Index df = (rank >= 0) ? n - X.cols() : n - rank;
    double s = resid.norm() / std::sqrt(double(df));
    Eigen::VectorXd se = s * ans.se();
    standard_error = se(1);
    double rss = resid.squaredNorm();
    double mss = (ans.fitted().array() - ans.fitted().mean()).pow(2).sum();
    r2 = mss / (mss + rss);
    long df_int = intercept; // 0 false 1 true
    r2_adjust = 1.0 - (1.0 - r2) * (static_cast<double>(n - df_int) / df);
    double tval = coeff / standard_error;
    p_value = misc::calc_tprob(tval, n);
}

}
