/*
 * regression.h
 *
 *  Created on: 5 Sep 2016
 *      Author: shingwanchoi
 */

#ifndef PRSICE_REGRESSION_H_
#define PRSICE_REGRESSION_H_

#include <Eigen/Dense>
#include <math.h>
#include <Rmath.h>

namespace Regression{
	void lm(Eigen::VectorXd y, Eigen::MatrixXd x, double &p_value, double &r2, double &r2_adjust, size_t thread=1,  bool intercept=true);
}

#endif /* PRSICE_REGRESSION_H_ */
