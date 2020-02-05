#ifndef GLM_HPP
#define GLM_HPP
#include "dcdflib.h"
#include "family.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <stdexcept>
template <typename family>
class GLM
{
public:
    GLM(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y,
        const Eigen::VectorXd& weights, const family& fam, double tol = 1e-8,
        int maxit = 100)
        : m_X(X)
        , m_Y(Y)
        , m_weights(weights)
        , m_nvars(X.cols())
        , m_nobs(X.rows())
        , m_family(fam)
        , m_beta(X.cols())
        , m_beta_prev(X.cols())
        , m_eta(X.rows())
        , m_var_mu(X.rows())
        , m_mu_eta(X.rows())
        , m_mu(X.rows())
        , m_z(X.rows())
        , m_w(weights)
        , m_se(X.cols())
        , m_maxit(maxit)
        , m_tol(tol)
    {
    }
    GLM(const Eigen::MatrixXd& X, const Eigen::VectorXd& Y, const family& fam,
        double tol = 1e-8, int maxit = 100)
        : m_X(X)
        , m_Y(Y)
        , m_weights(Eigen::VectorXd::Constant(X.rows(), 1))
        , m_nvars(X.cols())
        , m_nobs(X.rows())
        , m_family(fam)
        , m_beta(X.cols())
        , m_beta_prev(X.cols())
        , m_eta(X.rows())
        , m_var_mu(X.rows())
        , m_mu_eta(X.rows())
        , m_mu(X.rows())
        , m_z(X.rows())
        , m_w(Eigen::VectorXd::Constant(X.rows(), 1))
        , m_se(X.cols())
        , m_tol(tol)
        , m_maxit(maxit)
    {
    }
    virtual ~GLM() {}


    void init_parms(int type)
    {
        m_type = type;
        m_beta = Eigen::VectorXd::Zero(m_X.cols());
        m_eta = m_family.link(m_family.initialize(m_Y, m_w));
        m_mu = m_family.linkinv(m_eta);
        if (!m_family.validmu(m_mu) && !m_family.valideta(m_eta))
        {
            throw std::runtime_error("Error: GLM cannot find valid starting "
                                     "values");
        }
        // m_offset = Eigen::VectorXd::Zero(m_Y.rows());
        update_dev_resids();
        m_rank = m_nvars;
    }
    void init_parms()
    {
        // while type=2 should in theory be faster in most situation, it seems
        // to struggle when there is rank deficient cases and I haven't figure
        // out how to solve that
        // The problem is that it will both give a wrong beta and se when
        // there's rank deficiency. Now force back to the traditional QR
        // decomoposition
        //        size_t qr_time = 2 * m_nobs * m_nvars * m_nvars
        //                         - (2.0 / 3.0) * m_nvars * m_nvars * m_nvars;
        //        size_t fastqr_time = m_nobs * m_nvars * m_nvars
        //                             + (4.0 / 3.0) * m_nvars * m_nvars *
        //                             m_nvars;
        //        if (qr_time <= fastqr_time) { m_type = 1; }
        //        else
        //        {
        //            m_type = 2;
        //        }
        m_type = 1;
        m_beta = Eigen::VectorXd::Zero(m_X.cols());
        m_eta = m_family.link(m_family.initialize(m_Y, m_w));
        m_mu = m_family.linkinv(m_eta);
        if (!m_family.validmu(m_mu) && !m_family.valideta(m_eta))
        {
            throw std::runtime_error("Error: GLM cannot find valid starting "
                                     "values");
        }
        // m_offset = Eigen::VectorXd::Zero(m_Y.rows());
        update_dev_resids();
        m_rank = m_nvars;
    }
    int solve(int maxit = 100)
    {
        int i = 0;
        for (; i < maxit; ++i)
        {
            update_var_mu();
            update_mu_eta();
            update_z();
            update_w();
            solve_wls();
            update_eta();
            update_mu();
            update_dev_resids();
            run_step_halving(i);
            if (std::isinf(m_dev) && i == 0)
            {
                throw std::runtime_error("Error: cannot find valid starting "
                                         "values: please specify some");
            }
            if (converged())
            {
                m_converged = true;
                break;
            }
        }
        save_se();
        return std::min(i + 1, maxit);
    }
    const Eigen::VectorXd& get_beta() const { return m_beta; }
    const Eigen::VectorXd& get_se() const { return m_se; }
    double deviance() const { return m_dev; }
    bool has_converged() const { return m_converged; }
    double get_r2() const
    {
        double nulldev = m_family.dev_resids_sum(
            m_Y, Eigen::VectorXd::Constant(m_nobs, m_Y.sum() / m_nobs),
            m_weights);
        return (1.0 - std::exp((m_dev - nulldev) / static_cast<double>(m_nobs)))
               / (1.0 - std::exp(-nulldev / static_cast<double>(m_nobs)));
    }
    void get_stat(Eigen::Index idx, double& p, double& coeff, double& se) const
    {
        coeff = m_beta(idx);
        se = m_se(idx);
        double tvalue = coeff / se;
        p = chiprob_p(tvalue * tvalue, 1);
    }

private:
    const Eigen::MatrixXd m_X;
    const Eigen::VectorXd m_Y;
    const Eigen::VectorXd m_weights;
    const Eigen::Index m_nvars;
    const Eigen::Index m_nobs;
    family m_family;
    Eigen::LLT<Eigen::MatrixXd> m_Ch;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> m_PQR;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType m_Pmat;
    Eigen::MatrixXd m_Rinv;
    Eigen::VectorXd m_beta;
    Eigen::VectorXd m_beta_prev;
    Eigen::VectorXd m_eta;
    Eigen::VectorXd m_var_mu;
    Eigen::VectorXd m_mu_eta;
    Eigen::VectorXd m_mu;
    Eigen::VectorXd m_z;
    Eigen::VectorXd m_w;
    Eigen::VectorXd m_se;
    Eigen::VectorXd m_effects;
    // Eigen::VectorXd m_offset;
    double m_dev, m_devold;
    double m_tol = 1e-8;
    Eigen::Index m_rank;
    int m_maxit = 100;
    int m_type = 2;
    bool m_converged = false;

    bool converged() const
    {
        return (std::fabs(m_dev - m_devold) / (0.1 + std::fabs(m_dev)) < m_tol);
    }
    void update_eta() { m_eta = m_X * m_beta; }
    void update_var_mu() { m_var_mu = m_family.variance(m_mu); }
    void update_mu_eta() { m_mu_eta = m_family.mu_eta(m_eta); }
    void update_mu() { m_mu = m_family.linkinv(m_eta); }
    void update_z()
    {
        //   m_z = (m_eta.array() - m_offset.array())
        //       + (m_Y - m_mu).array() / m_mu_eta.array();
        m_z = m_eta.array() + (m_Y - m_mu).array() / m_mu_eta.array();
    }
    void update_w()
    {
        m_w = (m_weights.array() * m_mu_eta.array().square() / m_var_mu.array())
                  .array()
                  .sqrt();
    }
    void step_halve()
    {
        m_beta = 0.5 * (m_beta.array() + m_beta_prev.array());
        update_eta();
        update_mu();
    }
    void run_step_halving(int& iterr)
    {
        // check for infinite deviance
        if (std::isinf(m_dev))
        {
            int itrr = 0;
            while (std::isinf(m_dev))
            {
                ++itrr;
                if (itrr > m_maxit) break;
                step_halve();
                update_dev_resids_no_update();
            }
        }
        // check for boundary violations
        if (!m_family.valideta(m_eta) && m_family.validmu(m_mu))
        {
            int itrr = 0;
            while (!m_family.valideta(m_eta) && m_family.validmu(m_mu))
            {
                ++itrr;
                if (itrr > m_maxit) break;
                step_halve();
            }
            update_dev_resids_no_update();
        }
        if ((m_dev - m_devold) / (0.1 + std::abs(m_dev)) >= m_tol && iterr > 0)
        {
            int itrr = 0;
            while ((m_dev - m_devold) / (0.1 + std::abs(m_dev)) >= -m_tol)
            {
                ++itrr;
                if (itrr > m_maxit) break;
                step_halve();
                update_dev_resids_no_update();
            }
        }
    }
    void update_dev_resids_no_update()
    {
        m_dev = m_family.dev_resids_sum(m_Y, m_mu, m_weights);
    }
    void update_dev_resids()
    {
        m_devold = m_dev;
        m_dev = m_family.dev_resids_sum(m_Y, m_mu, m_weights);
    }

    Eigen::MatrixXd XtWX() const
    {
        return Eigen::MatrixXd(m_nvars, m_nvars)
            .setZero()
            .selfadjointView<Eigen::Lower>()
            .rankUpdate((m_w.asDiagonal() * m_X).adjoint());
    }
    void solve_wls()
    {
        m_beta_prev = m_beta;
        if (m_type == 0)
        {
            // use LLT
            m_Ch.compute(static_cast<Eigen::MatrixXd>(XtWX())
                             .selfadjointView<Eigen::Lower>());
            m_beta = m_Ch.solve((m_w.asDiagonal() * m_X).adjoint()
                                * (m_z.array() * m_w.array()).matrix());
        }
        else if (m_type == 1)
        {
            // use Col QR
            m_PQR.compute(m_w.asDiagonal() * m_X); // decompose the model matrix
            m_Pmat = (m_PQR.colsPermutation());
            m_rank = m_PQR.rank();
            if (m_rank == m_nvars)
            { // full rank case
                m_beta = m_PQR.solve((m_z.array() * m_w.array()).matrix());
            }
            else
            {
                m_Rinv = static_cast<Eigen::MatrixXd>(
                             m_PQR.matrixQR().topLeftCorner(m_rank, m_rank))
                             .triangularView<Eigen::Upper>()
                             .solve(Eigen::MatrixXd::Identity(m_rank, m_rank));
                m_effects = m_PQR.householderQ().adjoint()
                            * (m_z.array() * m_w.array()).matrix();
                m_beta.head(m_rank) = m_Rinv * m_effects.head(m_rank);
                m_beta = m_Pmat * m_beta;
                // create fitted values from effects
                // (can't use X*m_coef if X is rank-deficient)
                m_effects.tail(m_nobs - m_rank).setZero();
            }
        }
        else if (m_type == 2)
        {
            m_PQR.compute(static_cast<Eigen::MatrixXd>(XtWX())
                              .selfadjointView<Eigen::Lower>());
            m_Pmat = m_PQR.colsPermutation();
            m_rank = m_PQR.rank();
            if (m_rank == m_nvars)
            { // full rank case
                m_beta = m_PQR.solve((m_w.asDiagonal() * m_X).adjoint()
                                     * (m_z.array() * m_w.array()).matrix());
            }
            else
            {
                m_Rinv =
                    (Eigen::MatrixXd(
                         m_PQR.matrixQR().topLeftCorner(m_rank, m_rank))
                         .triangularView<Eigen::Upper>()
                         .solve(Eigen::MatrixXd::Identity(m_rank, m_rank)));
                m_effects = m_PQR.householderQ().adjoint()
                            * (m_w.asDiagonal() * m_X).adjoint()
                            * (m_z.array() * m_w.array()).matrix();
                m_beta.head(m_rank) = m_Rinv * m_effects.head(m_rank);
                m_beta = m_Pmat * m_beta;
                // create fitted values from effects
                // (can't use X*m_coef if X is rank-deficient)
                m_effects.tail(m_nvars - m_rank).setZero();
            }
        }
    }
    void save_se()
    {
        if (m_type == 0)
        {
            m_se = m_Ch.matrixL()
                       .solve(Eigen::MatrixXd::Identity(m_nvars, m_nvars))
                       .colwise()
                       .norm();
        }
        else if (m_type == 1 || m_type == 2)
        {
            if (m_rank == m_nvars)
            { // full rank case
                m_se = m_Pmat
                       * static_cast<Eigen::MatrixXd>(
                             m_PQR.matrixQR().topRows(m_nvars))
                             .triangularView<Eigen::Upper>()
                             .solve(Eigen::MatrixXd::Identity(m_nvars, m_nvars))
                             .rowwise()
                             .norm();
            }
            else
            {
                m_se.head(m_rank) = m_Rinv.rowwise().norm();
                m_se = m_Pmat * m_se;
            }
        }
        if (m_type == 2) { m_se = m_se.array().sqrt(); }
    }
};

#endif // GLM_NEW_HPP
