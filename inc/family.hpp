#ifndef FAMILY_HPP
#define FAMILY_HPP
#include <Eigen/Dense>
class Family
{
public:
    Family() {}
    virtual ~Family() {}

protected:
    double y_log_y(const double y, const double mu) const
    {
        return (y != 0.) ? (y * log(y / mu)) : 0;
    }
};

class Binomial : public Family
{
public:
    Binomial() {}
    virtual ~Binomial() {}
    Eigen::VectorXd variance(const Eigen::VectorXd& eta) const
    {
        return eta.array() * (1 - eta.array());
    }
    Eigen::VectorXd mu_eta(const Eigen::VectorXd& eta) const
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
    Eigen::VectorXd linkinv(const Eigen::VectorXd& eta) const
    {
        Eigen::VectorXd ans = eta;
        long n = eta.rows();
        double etai, temp;
        for (long i = 0; i < n; ++i)
        {
            etai = eta(i);
            temp =
                (etai < -30)
                    ? std::numeric_limits<double>::epsilon()
                    : ((etai > 30) ? 1 / std::numeric_limits<double>::epsilon()
                                   : exp(etai));
            ans(i) = temp / (1.0 + temp);
        }
        return ans;
    }
    Eigen::VectorXd initialize(const Eigen::VectorXd& y,
                               const Eigen::VectorXd weights) const
    {
        // need to add script to "zero out" y when weight == 0
        return ((weights.array() * y.array()).array() + 0.5)
               / (weights.array() + 1);
    }
    Eigen::VectorXd link(const Eigen::VectorXd& mu_start) const
    {
        return mu_start.array().log() - (1 - mu_start.array()).log();
    }

    double dev_resids_sum(const Eigen::VectorXd& y, const Eigen::VectorXd& mu,
                          const Eigen::VectorXd& wt) const
    {
        // wt should always be 1 in our case, because we don't allow weighted
        // regression
        const Eigen::Index n = y.rows();
        const Eigen::Index lmu = mu.rows();
        const Eigen::Index lwt = wt.rows();
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
                // do the sum here directly
                ans += 2 * wt[(lwt > 1) ? i : 0]
                       * (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
            }
        }
        else
        {
            mui = mu[0];
            for (Eigen::Index i = 0; i < n; ++i)
            {
                yi = y(i);
                ans += 2 * wt[(lwt > 1) ? i : 0]
                       * (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
            }
        }
        return ans;
    }
    double var_res() { return 1; }
    bool valideta(const Eigen::VectorXd& /*eta*/) const { return true; }
    bool validmu(const Eigen::VectorXd& mu) const
    {
        return (mu.array().isFinite() == 1).all()
               && (mu.array() > 0.0 && mu.array() < 1.0).all();
    }
};

#endif // FAMILY_HPP
