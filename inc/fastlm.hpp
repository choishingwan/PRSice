#ifndef FASTLM_HPP
#define FASTLM_HPP
#include <Eigen/Dense>
#include <limits>
// From RcppEigen
class lm
{
public:
    lm(const Eigen::MatrixXd&, const Eigen::VectorXd&);
    lm() {}
    inline Eigen::ArrayXd Dplus(const Eigen::ArrayXd& d)
    {
        Eigen::ArrayXd di(d.size());
        double comp(d.maxCoeff() * threshold());
        for (int j = 0; j < d.size(); ++j)
            di[j] = (d[j] < comp) ? 0. : 1. / d[j];
        m_r = (di != 0.).count();
        return di;
    }
    Eigen::MatrixXd I_p() const { return Eigen::MatrixXd::Identity(m_p, m_p); }
    Eigen::MatrixXd XtX() const
    {
        return Eigen::MatrixXd(m_p, m_p)
            .setZero()
            .selfadjointView<Eigen::Lower>()
            .rankUpdate(m_X.adjoint());
    }
    Eigen::MatrixXd::RealScalar threshold() const
    {
        return m_usePrescribedThreshold
                   ? m_prescribedThreshold
                   : std::numeric_limits<double>::epsilon() * m_p;
    }
    const Eigen::VectorXd& se() const { return m_se; }
    const Eigen::VectorXd& coef() const { return m_coef; }
    const Eigen::VectorXd& fitted() const { return m_fitted; }
    int rank() const { return m_r; }
    // return a lm object without copying
    lm& setThreshold(const Eigen::MatrixXd::RealScalar&);

protected:
    Eigen::MatrixXd m_X;
    Eigen::VectorXd m_y;
    Eigen::Index m_n;
    Eigen::Index m_p;
    Eigen::VectorXd m_coef;
    int m_r;
    Eigen::VectorXd m_fitted;
    Eigen::VectorXd m_se;
    Eigen::MatrixXd::RealScalar m_prescribedThreshold;
    bool m_usePrescribedThreshold;
};

class ColPivQR : public lm
{
public:
    ColPivQR(const Eigen::MatrixXd&, const Eigen::VectorXd&);
};

class Llt : public lm
{
public:
    Llt(const Eigen::MatrixXd&, const Eigen::VectorXd&);
};

class Ldlt : public lm
{
public:
    Ldlt(const Eigen::MatrixXd&, const Eigen::VectorXd&);
};

class QR : public lm
{
public:
    QR(const Eigen::MatrixXd&, const Eigen::VectorXd&);
};

class SVD : public lm
{
public:
    SVD(const Eigen::MatrixXd&, const Eigen::VectorXd&);
};

class SymmEigen : public lm
{
public:
    SymmEigen(const Eigen::MatrixXd&, const Eigen::VectorXd&);
};


#endif // FASTLM_HPP
