/*
 * regression.cpp
 *
 *  Created on: 5 Sep 2016
 *      Author: shingwanchoi
 */

#include "regression.h"
namespace Regression{
	void lm(Eigen::VectorXd &y, Eigen::MatrixXd &x, double &p_value, double &r2, double &r2_adjust, size_t thread, bool intercept){
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
			rss+=residual*residual;
		}
		int rank = z.rank();
		int n = x.cols();
		int rdf = n-rank;
		double resvar = rss/(double)rdf;
		int df_int = intercept; //0 false 1 true
		r2 = mss/(mss+rss);
		r2_adjust = 1.0- (1.0-r2)/((double)(n - df_int)/(double)rdf);
		double fstat = mss/(double)(rank-df_int);
		double df = rank-df_int;
		double dendf = rdf;
		p_value=pf(fstat, df,dendf,false);
	}
}



