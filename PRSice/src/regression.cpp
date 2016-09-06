/*
 * regression.cpp
 *
 *  Created on: 5 Sep 2016
 *      Author: shingwanchoi
 */

#include "regression.h"
namespace Regression{

	void linear_regression(Eigen::VectorXd &y, Eigen::MatrixXd &x, double &p_value, double &r2, double &r2_adjust, size_t thread, bool intercept){
		if(intercept){
			x.conservativeResize(x.rows(), x.cols()+1);
			x.col(x.cols()-1) = Eigen::VectorXd::Constant(x.rows(),1.0);
		}
		Eigen::ColPivHouseholderQR<Eigen::MatrixXd> z(x);
		Eigen::VectorXd beta = z.solve(y);
		Eigen::MatrixXd fitted = x*beta;
		Eigen::VectorXd residual = y-fitted;
		double mss = 0.0;
		double rss = 0.0;
		for(size_t i = 0; i < x.rows(); ++i){
			mss+=pow(fitted(i)-fitted.mean(),2);
			rss+=residual(i)*residual(i);
		}
		int rank = z.rank();
		int n = x.rows();
		int rdf = n-rank;
		double resvar = rss/(double)rdf;
		int df_int = intercept; //0 false 1 true
		r2 = mss/(mss+rss);
		r2_adjust = 1.0- (1.0-r2)*((double)(n - df_int)/(double)rdf);
		Eigen::MatrixXd R = z.matrixR().topLeftCorner(rank, rank).triangularView<Eigen::Upper>();
		Eigen::VectorXd se = ((R.transpose()*R).inverse().diagonal()*resvar).array().sqrt();
		// Remember, only the coefficient's order is wrong e.g. intercept at the end
		Eigen::VectorXd est = beta;
		double tval = beta(0)/se(intercept); // only interested in the one coefficient
		boost::math::students_t dist(rdf);
		p_value = 2*boost::math::cdf(boost::math::complement(dist, fabs(tval)));
	}
}



