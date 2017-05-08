/*
 * regression.h
 *
 *  Created on: 5 Sep 2016
 *      Author: shingwanchoi
 */

#ifndef PRSICE_REGRESSION_H_
#define PRSICE_REGRESSION_H_

#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/normal.hpp>
#include <cstdio>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <stdexcept>

namespace Regression{
	void linear_regression(const Eigen::VectorXd &y, const Eigen::MatrixXd &A, double &p_value, double &r2, double &r2_adjust, double &coeff, size_t thread=1, bool intercept=true);
	void glm(const Eigen::VectorXd &y, const Eigen::MatrixXd &x, double &p_value, double &r2 , double &coeff, size_t max_iter=25,  size_t thread=1, bool intercept=true);
	Eigen::VectorXd logit_linkinv(const Eigen::VectorXd &eta);
	Eigen::VectorXd logit_variance(const Eigen::VectorXd &eta);
	Eigen::VectorXd logit_mu_eta(const Eigen::VectorXd &eta);
	Eigen::VectorXd binomial_dev_resids(const Eigen::VectorXd &y, const Eigen::VectorXd &mu, const Eigen::VectorXd &wt);
	double binomial_dev_resids_sum(const Eigen::VectorXd &y, const Eigen::VectorXd &mu, const Eigen::VectorXd &wt);
	void logit_both(const Eigen::VectorXd &eta, Eigen::VectorXd &g, Eigen::VectorXd &gprime);
	inline double y_log_y(double y, double mu){
	    return (y != 0.) ? (y * log(y/mu)) : 0;
	}
}

#endif /* PRSICE_REGRESSION_H_ */
