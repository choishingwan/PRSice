#include "fastlm.hpp"


lm::lm(const Eigen::MatrixXd& X, const Eigen::VectorXd& y)
    : m_X(X)
    , m_y(y)
    , m_n(X.rows())
    , m_p(X.cols())
    , m_coef(Eigen::VectorXd::Constant(
          m_p, std::numeric_limits<double>::quiet_NaN()))
    , m_r(-1)
    , m_fitted(m_n)
    , m_se(Eigen::VectorXd::Constant(m_p,
                                     std::numeric_limits<double>::quiet_NaN()))
    , m_usePrescribedThreshold(false)
{
}

lm& lm::setThreshold(const Eigen::MatrixXd::RealScalar& threshold)
{
    m_usePrescribedThreshold = true;
    m_prescribedThreshold = threshold;
    return *this;
}
ColPivQR::ColPivQR(const Eigen::MatrixXd& X, const Eigen::VectorXd& y)
    : lm(X, y)
{
    // decompose the model matrix
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> PQR(X);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType Pmat(
        PQR.colsPermutation());
    m_r = PQR.rank();
    if (m_r == m_p)
    { // full rank case
        m_coef = PQR.solve(y);
        m_fitted = X * m_coef;
        m_se = Pmat
               * PQR.matrixQR()
                     .topRows(m_p)
                     .triangularView<Eigen::Upper>()
                     .solve(I_p())
                     .rowwise()
                     .norm();
        return;
    }
    Eigen::MatrixXd Rinv(PQR.matrixQR()
                             .topLeftCorner(m_r, m_r)
                             .triangularView<Eigen::Upper>()
                             .solve(Eigen::MatrixXd::Identity(m_r, m_r)));
    m_se.head(m_r) = Rinv.rowwise().norm();
    m_se = Pmat * m_se;

    Eigen::VectorXd effects(PQR.householderQ().adjoint() * y);
    m_coef.head(m_r) = Rinv * effects.head(m_r);
    m_coef = Pmat * m_coef;
    // create fitted values from effects
    // (can't use X*m_coef if X is rank-deficient)
    effects.tail(m_n - m_r).setZero();
    m_fitted = PQR.householderQ() * effects;
}

QR::QR(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) : lm(X, y)
{
    Eigen::HouseholderQR<Eigen::MatrixXd> QR(X);
    m_coef = QR.solve(y);
    m_fitted = X * m_coef;
    m_se = QR.matrixQR()
               .topRows(m_p)
               .triangularView<Eigen::Upper>()
               .solve(I_p())
               .rowwise()
               .norm();
}


Llt::Llt(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) : lm(X, y)
{
    Eigen::LLT<Eigen::MatrixXd> Ch(XtX().selfadjointView<Eigen::Lower>());
    m_coef = Ch.solve(X.adjoint() * y);
    m_fitted = X * m_coef;
    m_se = Ch.matrixL().solve(I_p()).colwise().norm();
}

Ldlt::Ldlt(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) : lm(X, y)
{
    Eigen::LDLT<Eigen::MatrixXd> Ch(XtX().selfadjointView<Eigen::Lower>());
    Dplus(Ch.vectorD()); // to set the rank
    // FIXME: Check on the permutation in the LDLT and incorporate it in
    // the coefficients and the standard error computation.
    //	m_coef            = Ch.matrixL().adjoint().
    //	    solve(Dplus(D) * Ch.matrixL().solve(X.adjoint() * y));
    m_coef = Ch.solve(X.adjoint() * y);
    m_fitted = X * m_coef;
    m_se = Ch.solve(I_p()).diagonal().array().sqrt();
}


SVD::SVD(const Eigen::MatrixXd& X, const Eigen::VectorXd& y) : lm(X, y)
{
    Eigen::JacobiSVD<Eigen::MatrixXd> UDV(
        X.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV));
    Eigen::MatrixXd VDi(
        UDV.matrixV()
        * Dplus(UDV.singularValues().array()).matrix().asDiagonal());
    m_coef = VDi * UDV.matrixU().adjoint() * y;
    m_fitted = X * m_coef;
    m_se = VDi.rowwise().norm();
}

SymmEigen::SymmEigen(const Eigen::MatrixXd& X, const Eigen::VectorXd& y)
    : lm(X, y)
{
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(
        XtX().selfadjointView<Eigen::Lower>());
    Eigen::MatrixXd VDi(
        eig.eigenvectors()
        * Dplus(eig.eigenvalues().array()).sqrt().matrix().asDiagonal());
    m_coef = VDi * VDi.adjoint() * X.adjoint() * y;
    m_fitted = X * m_coef;
    m_se = VDi.rowwise().norm();
}
