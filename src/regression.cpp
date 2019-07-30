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

    long n = A.rows();
    long rank = z.rank();
    long rdf = n - rank;
    long df_int = intercept; // 0 false 1 true
    r2_adjust =
        1.0
        - (1.0 - r2)
              * (static_cast<double>(n - df_int) / static_cast<double>(rdf));


    double resvar = rss / static_cast<double>(rdf);
    long se_index = intercept;
    for (long ind = 0; ind < beta.rows(); ++ind)
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

Eigen::VectorXd logit_variance(const Eigen::VectorXd& eta)
{
    // Eigen::VectorXd ans = eta;
    return eta.array() * (1 - eta.array());
    // for (long i = 0; i < eta.rows(); ++i) ans(i) = eta(i) * (1.0 - eta(i));
    // return ans;
}

Eigen::VectorXd logit_mu_eta(const Eigen::VectorXd& eta)
{
    Eigen::VectorXd ans = eta;
    long n = eta.rows();
    double etai, opexp;
    const double limit = std::numeric_limits<double>::epsilon();
    for (long i = 0; i < n; ++i)
    {
        etai = eta(i);
        opexp = 1 + exp(etai);
        ans(i) =
            (etai > 30 || etai < -30) ? limit : exp(etai) / (opexp * opexp);
    }
    return ans;
}

Eigen::VectorXd logit_linkinv(const Eigen::VectorXd& eta)
{
    Eigen::VectorXd ans = eta;
    long n = eta.rows();
    double etai, temp;
    for (long i = 0; i < n; ++i)
    {
        etai = eta(i);
        temp = (etai < -30)
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
    long n = eta.rows();
    g = eta;
    gprime = eta;
    double etai, temp, opexp;
    const double limit = std::numeric_limits<double>::epsilon();
    for (long i = 0; i < n; ++i)
    {
        etai = eta(i);
        temp = (etai < -30) ? limit : ((etai > 30) ? 1 / limit : exp(etai));
        g(i) = temp / (1.0 + temp);
        opexp = 1 + exp(etai);
        gprime(i) =
            (etai > 30 || etai < -30) ? limit : exp(etai) / (opexp * opexp);
    }
}

Eigen::VectorXd binomial_dev_resids(const Eigen::VectorXd& y,
                                    const Eigen::VectorXd& mu,
                                    const Eigen::VectorXd& wt)
{
    Eigen::Index n = y.rows();
    Eigen::Index lmu = mu.rows(), lwt = wt.rows();
    Eigen::VectorXd ans = y;
    if (lmu != n && lmu != 1)
    {
        std::string error_message =
            "Argument mu must be a numeric vector of length 1 or length "
            + std::to_string(n);
        throw std::runtime_error(error_message);
    }
    if (lwt != n && lwt != 1)
    {
        std::string error_message =
            "Argument wt must be a numeric vector of length 1 or length "
            + std::to_string(n);
        throw std::runtime_error(error_message);
    }
    double mui, yi;
    if (lmu > 1)
    {
        for (Eigen::Index i = 0; i < n; ++i)
        {
            mui = mu(i);
            yi = y(i);
            ans(i) = 2 * wt((lwt > 1) ? i : 0)
                     * (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
        }
    }
    else
    {
        mui = mu[0];
        for (Eigen::Index i = 0; i < n; ++i)
        {
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
    // wt should always be 1 in our case, because we don't allow weighted
    // regression
    const Eigen::Index n = y.rows();
    const Eigen::Index lmu = mu.rows(), lwt = wt.rows();
    double ans = 0.0;
    if (lmu != n && lmu != 1)
    {
        std::string error_message =
            "Argument mu must be a numeric vector of length 1 or length "
            + std::to_string(n);
        throw std::runtime_error(error_message);
    }
    if (lwt != n && lwt != 1)
    {
        std::string error_message =
            "Argument wt must be a numeric vector of length 1 or length "
            + std::to_string(n);
        throw std::runtime_error(error_message);
    }
    double mui, yi;
    if (lmu > 1)
    {
        for (Eigen::Index i = 0; i < n; ++i)
        {
            mui = mu(i);
            yi = y(i);
            ans += 2 * (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
            // ans += 2 * wt((lwt > 1) ? i : 0) * (y_log_y(yi, mui) + y_log_y(1
            // - yi, 1 - mui));
        }
    }
    else
    {
        mui = mu[0];
        for (Eigen::Index i = 0; i < n; ++i)
        {
            yi = y(i);
            ans += 2 * (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
            // ans += 2 * wt((lwt > 1) ? i : 0) * (y_log_y(yi, mui) + y_log_y(1
            // - yi, 1 - mui));
        }
    }
    return ans;
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
    /*
    const Eigen::Index nobs = y.rows();
    const Eigen::Index nvars = x.cols();
    Eigen::VectorXd weights = Eigen::VectorXd::Ones(nobs);
    Eigen::VectorXd eta = family.link(family.initialize(y, weights));
    Eigen::VectorXd mu = family.linkinv(eta);
    Eigen::MatrixXd A = x;
    double devold = family.dev_resids_sum(y, mu, weights);
    double dev = 0.0;
    // Iterative reweighting
    Eigen::VectorXd varmu;
    Eigen::VectorXd mu_eta_val;
    Eigen::VectorXd z = weights;
    Eigen::VectorXd w = weights;
    Eigen::VectorXd beta;
    bool good;
    Eigen::VectorXd fit;
    Eigen::VectorXd start;
    bool converge = false;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr;
    qr.setThreshold(
        std::min(1e-7, std::numeric_limits<double>::epsilon() / 1000));
    double mu_eta_val_store, mu_store;
    for (size_t iter = 0; iter < max_iter; ++iter)
    {
        mu_eta_val = family.mu_eta(eta);
        // we can ignore weight array check because we don't allow weighted glm
        // at the moment. So all weight must be 1
        //			good = (weights.array()>0 && mu_eta_val.array() !=
        // 0).select(weights.array(), 0);
        z = eta.array() + (y - mu).array() / mu_eta_val.array();
        w = (weights.array() * mu_eta_val.array().square()
             / family.variance(mu).array())
                .array()
                .sqrt();
        qr.compute(w.asDiagonal() * x);
        beta = qr.solve((z.array() * w.array()).matrix());
        if (nobs < qr.rank())
        {
            std::string error_message =
                "X matrix has rank " + std::to_string(qr.rank()) + "but only "
                + std::to_string(nobs) + " observations";
            throw std::runtime_error(error_message);
        }
        eta = A * start;
        mu = family.linkinv(eta);
        dev = family.dev_resids_sum(y, mu, weights);
        // R only use 1e-8 here
        if (fabs(dev - devold) / (0.1 + fabs(dev)) < 1e-8)
        {
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
    Eigen::VectorXd wtdmu =
        (intercept) ? Eigen::VectorXd::Constant(nobs, y.sum() / nobs)
                    : family.linkinv(Eigen::VectorXd::Zero(nobs));
    double nulldev = family.dev_resids_sum(y, wtdmu, weights);
    Eigen::Index rank = qr.rank();
    Eigen::MatrixXd R =
        qr.matrixQR().topLeftCorner(rank, rank).triangularView<Eigen::Upper>();
    Eigen::VectorXd se =
        ((R.transpose() * R).inverse().diagonal()).array().sqrt();
    // again, we are ony interested in one of the variable
    r2 = (1.0 - std::exp((dev - nulldev) / static_cast<double>(nobs)))
         / (1 - std::exp(-nulldev / static_cast<double>(nobs)));
    Eigen::Index se_index = intercept;
    for (Eigen::Index ind = 0; ind < start.rows(); ++ind)
    {
        if (qr.colsPermutation().indices()(ind) == intercept)
        {
            se_index = ind;
            break;
        }
    }

    double tvalue = start(intercept) / se(se_index);
    coeff = start(intercept);
    p_value = chiprob_p(tvalue * tvalue, 1);
    standard_error = se(se_index);
    */
}
}
