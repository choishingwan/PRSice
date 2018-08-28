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


#include "regression.hpp"
namespace Regression
{

// on purposely perform the copying of x
void linear_regression(const Eigen::VectorXd& y, const Eigen::MatrixXd& A,
                       double& p_value, double& r2, double& r2_adjust,
                       double& coeff, double& standard_error, intptr_t thread,
                       bool intercept)
{
    Eigen::setNbThreads(thread);
    // in more general cases, the following is needed (adding the intercept)
    // but here in PRSice, we always include the intercept, so we will skip
    // this to speed things up
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> z(A);
    Eigen::VectorXd beta = z.solve(y);
    Eigen::VectorXd fitted = A * beta;
    double rss = (A * beta - y).squaredNorm();
    double mss = (fitted.array() - fitted.mean()).pow(2).sum();
    r2 = mss / (mss + rss);

    int n = A.rows();
    int rank = z.rank();
    int rdf = n - rank;
    int df_int = intercept; // 0 false 1 true
    r2_adjust = 1.0 - (1.0 - r2) * ((double) (n - df_int) / (double) rdf);


    double resvar = rss / (double) rdf;
    size_t se_index = intercept;
    for (size_t ind = 0; ind < beta.rows(); ++ind) {
        if (z.colsPermutation().indices()(ind) == intercept) {
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

Eigen::VectorXd logit_variance(const Eigen::VectorXd& eta)
{
    Eigen::VectorXd ans = eta;
    for (size_t i = 0; i < eta.rows(); ++i) ans(i) = eta(i) * (1.0 - eta(i));
    return ans;
}

Eigen::VectorXd logit_mu_eta(const Eigen::VectorXd& eta)
{
    Eigen::VectorXd ans = eta;
    int n = eta.rows();
    for (size_t i = 0; i < n; ++i) {
        double etai = eta(i);
        double opexp = 1 + exp(etai);
        ans(i) = (etai > 30 || etai < -30)
                     ? std::numeric_limits<double>::epsilon()
                     : exp(etai) / (opexp * opexp);
    }
    return ans;
}

Eigen::VectorXd logit_linkinv(const Eigen::VectorXd& eta)
{
    Eigen::VectorXd ans = eta;
    int n = eta.rows();
    for (size_t i = 0; i < n; ++i) {
        double etai = eta(i);
        double temp =
            (etai < -30)
                ? std::numeric_limits<double>::epsilon()
                : ((etai > 30) ? 1 / std::numeric_limits<double>::epsilon()
                               : exp(etai));
        ans(i) = temp / (1.0 + temp);
    }
    return ans;
}

void logit_both(const Eigen::VectorXd& eta, Eigen::VectorXd& g,
                Eigen::VectorXd& gprime)
{
    int n = eta.rows();
    g = eta;
    gprime = eta;
    for (size_t i = 0; i < n; ++i) {
        double etai = eta(i);
        double temp =
            (etai < -30)
                ? std::numeric_limits<double>::epsilon()
                : ((etai > 30) ? 1 / std::numeric_limits<double>::epsilon()
                               : exp(etai));
        g(i) = temp / (1.0 + temp);
        double opexp = 1 + exp(etai);
        gprime(i) = (etai > 30 || etai < -30)
                        ? std::numeric_limits<double>::epsilon()
                        : exp(etai) / (opexp * opexp);
    }
}

Eigen::VectorXd binomial_dev_resids(const Eigen::VectorXd& y,
                                    const Eigen::VectorXd& mu,
                                    const Eigen::VectorXd& wt)
{
    Eigen::Index n = y.rows();
    Eigen::Index lmu = mu.rows(), lwt = wt.rows();
    Eigen::VectorXd ans = y;
    if (lmu != n && lmu != 1) {
        std::string error_message =
            "Argument mu must be a numeric vector of length 1 or length "
            + std::to_string(n);
        throw std::runtime_error(error_message);
    }
    if (lwt != n && lwt != 1) {
        std::string error_message =
            "Argument wt must be a numeric vector of length 1 or length "
            + std::to_string(n);
        throw std::runtime_error(error_message);
    }
    double mui, yi;
    if (lmu > 1) {
        for (Eigen::Index i = 0; i < n; ++i) {
            mui = mu(i);
            yi = y(i);
            ans(i) = 2 * wt((lwt > 1) ? i : 0)
                     * (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
        }
    }
    else
    {
        mui = mu[0];
        for (Eigen::Index i = 0; i < n; ++i) {
            yi = y(i);
            ans(i) = 2 * wt((lwt > 1) ? i : 0)
                     * (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
        }
    }
    return ans;
}

double binomial_dev_resids_sum(const Eigen::VectorXd& y,
                               const Eigen::VectorXd& mu,
                               const Eigen::VectorXd& wt)
{
    Eigen::Index n = y.rows();
    Eigen::Index lmu = mu.rows(), lwt = wt.rows();
    double ans = 0.0;
    if (lmu != n && lmu != 1) {
        std::string error_message =
            "Argument mu must be a numeric vector of length 1 or length "
            + std::to_string(n);
        throw std::runtime_error(error_message);
    }
    if (lwt != n && lwt != 1) {
        std::string error_message =
            "Argument wt must be a numeric vector of length 1 or length "
            + std::to_string(n);
        throw std::runtime_error(error_message);
    }
    double mui, yi;
    if (lmu > 1) {
        for (Eigen::Index i = 0; i < n; ++i) {
            mui = mu(i);
            yi = y(i);
            ans += 2 * wt((lwt > 1) ? i : 0)
                   * (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
        }
    }
    else
    {
        mui = mu[0];
        for (Eigen::Index i = 0; i < n; ++i) {
            yi = y(i);
            ans += 2 * wt((lwt > 1) ? i : 0)
                   * (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
        }
    }
    return ans;
}

// This is an unsafe version of R's glm.fit
// unsafe as in I have skipped some of the checking
void glm(const Eigen::VectorXd& y, const Eigen::MatrixXd& x, double& p_value,
         double& r2, double& coeff, double& standard_error, size_t max_iter,
         intptr_t thread, bool intercept)
{
    Eigen::setNbThreads(static_cast<int>(thread));
    Eigen::MatrixXd A = x;
    Eigen::Index nobs = y.rows();
    Eigen::Index nvars = A.cols();
    Eigen::VectorXd weights = Eigen::VectorXd::Ones(nobs);
    Eigen::VectorXd mustart =
        (weights.array() * y.array() + 0.5) / (weights.array() + 1);
    Eigen::VectorXd eta =
        (mustart.array() / (1 - mustart.array())).array().log();
    Eigen::VectorXd mu = logit_linkinv(eta);
    double devold = binomial_dev_resids_sum(y, mu, weights), dev = 0.0;
    // Iterative reweighting
    Eigen::MatrixXd A_tmp;
    Eigen::VectorXd varmu, mu_eta_val, z, w, good, fit, start;
    bool converge = false;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr;
    qr.setThreshold(
        std::min(1e-7, std::numeric_limits<double>::epsilon() / 1000));
    for (size_t iter = 0; iter < max_iter; ++iter) {
        // varmu = (weights.array()>0).select(mu.array()*(1-mu.array()),0);
        mu_eta_val = logit_mu_eta(eta);
        //			good = (weights.array()>0 && mu_eta_val.array() !=
        // 0).select(weights.array(), 0);
        good = (mu_eta_val.array() != 0).select(mu_eta_val.array(), 0);
        z = weights;
        w = weights;
        Eigen::Index i_good = 0;
        for (Eigen::Index i_weights = 0; i_weights < good.rows(); ++i_weights) {
            if (good(i_weights) > 0) {
                // because offset is 0, we ignore it
                z(i_good) =
                    eta(i_weights)
                    + (y(i_weights) - mu(i_weights)) / mu_eta_val(i_weights);
                w(i_good) =
                    std::sqrt(mu_eta_val(i_weights) * mu_eta_val(i_weights)
                              / (mu(i_weights) * (1 - mu(i_weights))));
                A.row(i_good) = A.row(i_weights);
                i_good++;
            }
        }
        z.conservativeResize(i_good);
        w.conservativeResize(i_good);
        A.conservativeResize(i_good, nvars);
        A_tmp = A;
        for (Eigen::Index i = 0; i < nvars; ++i) {
            A_tmp.col(i) = A.col(i).array() * w.array();
        }
        qr.compute(A_tmp);
        start = qr.solve(Eigen::MatrixXd(z.array() * w.array()));
        if (nobs < qr.rank()) {
            std::string error_message =
                "X matrix has rank " + std::to_string(qr.rank()) + "but only "
                + std::to_string(nobs) + " observations";
            throw std::runtime_error(error_message);
        }
        eta = A * start;
        mu = logit_linkinv(eta);
        dev = binomial_dev_resids_sum(y, mu, weights);
        // R only use 1e-8 here
        if (fabs(dev - devold) / (0.1 + fabs(dev)) < 1e-8) {
            converge = true;
            break;
        }
        else
        {
            devold = dev;
        }
    }
    r2 = 0;
    coeff = 0;
    p_value = -1;
    if (!converge) throw std::runtime_error("GLM algorithm did not converge");
    Eigen::VectorXd residuals =
        (y.array() - mu.array()) / (logit_mu_eta(eta).array());
    int sum_good = 0;
    int weight_bad = 0;
    for (Eigen::Index i = 0; i < good.rows(); ++i) {
        if (!misc::logically_equal(good(i), 0)) sum_good++;
        if (misc::logically_equal(weights(i), 0)) weight_bad++;
    }
    Eigen::VectorXd wtdmu =
        (intercept) ? Eigen::VectorXd::Constant(
                          nobs, (weights.array() * y.array()).array().sum()
                                    / weights.array().sum())
                    : logit_linkinv(Eigen::VectorXd::Zero(nobs));
    double nulldev = binomial_dev_resids_sum(y, wtdmu, weights);
    Eigen::Index rank = qr.rank();
    Eigen::MatrixXd R =
        qr.matrixQR().topLeftCorner(rank, rank).triangularView<Eigen::Upper>();
    Eigen::VectorXd se =
        ((R.transpose() * R).inverse().diagonal()).array().sqrt();
    // again, we are ony interested in one of the variable
    r2 = (1.0 - std::exp((dev - nulldev) / static_cast<double>(nobs)))
         / (1 - std::exp(-nulldev / static_cast<double>(nobs)));
    Eigen::Index se_index = intercept;
    for (Eigen::Index ind = 0; ind < start.rows(); ++ind) {
        if (qr.colsPermutation().indices()(ind) == intercept) {
            se_index = ind;
            break;
        }
    }

    double tvalue = start(intercept) / se(se_index);
    coeff = start(intercept);
    p_value = misc::chiprob_p(tvalue * tvalue, 1);
    standard_error = se(se_index);
}
}
