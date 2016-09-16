/*
 * regression.cpp
 *
 *  Created on: 5 Sep 2016
 *      Author: shingwanchoi
 */

#include "regression.h"
namespace Regression{

	// on purposely perform the copying of x
	void linear_regression(const Eigen::VectorXd &y, const Eigen::MatrixXd &A, double &p_value, double &r2, double &r2_adjust, size_t thread, bool intercept){
		Eigen::setNbThreads(thread);
		// in more general cases, the following is needed (adding the intercept)
		// but here in PRSice, we always include the intercept, so we will skip
		// this to speed things up
		/*
		Eigen::MatrixXd A;
		if(intercept){
			A=Eigen::MatrixXd::Zero(x.rows(), x.cols()+1);
			A.col(0) = Eigen::VectorXd::Constant(x.rows(),1.0);
			A.block(0,1,x.rows(), x.cols()) = x;
		}
		else A = x;
		*/
		Eigen::ColPivHouseholderQR<Eigen::MatrixXd> z(A);
		Eigen::VectorXd beta = z.solve(y);
		Eigen::MatrixXd fitted = A*beta;
		Eigen::VectorXd residual = y-fitted;
		double mss = 0.0;
		double rss = 0.0;
		double fitted_mean = fitted.mean();
		for(size_t i = 0; i < A.rows(); ++i){
			mss+=pow(fitted(i)-fitted_mean,2);
			rss+=residual(i)*residual(i);
		}
		int rank = z.rank();
		int n = A.rows();
		int rdf = n-rank;
		double resvar = rss/(double)rdf;
		int df_int = intercept; //0 false 1 true
		r2 = mss/(mss+rss);
		r2_adjust = 1.0- (1.0-r2)*((double)(n - df_int)/(double)rdf);
		Eigen::MatrixXd R = z.matrixR().topLeftCorner(rank, rank).triangularView<Eigen::Upper>();
		Eigen::VectorXd se = ((R.transpose()*R).inverse().diagonal()*resvar).array().sqrt();
		// Remember, only the coefficient's order is wrong e.g. intercept at the end
		//Eigen::VectorXd est = beta.transpose()*Eigen::MatrixXd(z.colsPermutation());
		double tval = beta(intercept)/se(intercept); // only interested in the one coefficient
		boost::math::students_t dist(rdf);
		p_value = 2*boost::math::cdf(boost::math::complement(dist, fabs(tval)));
	}

	Eigen::VectorXd logit_variance(const Eigen::VectorXd &eta){
		Eigen::VectorXd ans = eta;
		for(size_t i = 0; i < eta.rows(); ++i) ans(i) = eta(i)*(1.0-eta(i));
		return ans;
	}

	Eigen::VectorXd logit_mu_eta(const Eigen::VectorXd &eta){
		Eigen::VectorXd ans = eta;
		int n = eta.rows();
		for(size_t i = 0; i < n; ++i){
			double etai =eta(i);
			double opexp = 1+exp(etai);
			ans(i) = (etai>30 || etai <-30)? std::numeric_limits<double>::epsilon():exp(etai)/(opexp*opexp);
		}
		return ans;
	}

	Eigen::VectorXd logit_linkinv(const Eigen::VectorXd &eta){
		Eigen::VectorXd ans = eta;
		int n = eta.rows();
		for(size_t i = 0; i < n; ++i){
			double etai = eta(i);
			double temp = (etai<-30)? std::numeric_limits<double>::epsilon():
					((etai>30)? 1/ std::numeric_limits<double>::epsilon() :exp(etai));
			ans(i) = temp/(1.0+temp);
		}
		return ans;
	}

	void logit_both(const Eigen::VectorXd &eta, Eigen::VectorXd &g, Eigen::VectorXd &gprime){
		int n = eta.rows();
		g=eta;
		gprime=eta;
		for(size_t i = 0; i < n; ++i){
			double etai = eta(i);
			double temp = (etai<-30)? std::numeric_limits<double>::epsilon():
					((etai>30)? 1/ std::numeric_limits<double>::epsilon() :exp(etai));
			g(i) = temp/(1.0+temp);
			double opexp = 1+exp(etai);
			gprime(i) = (etai>30 || etai <-30)? std::numeric_limits<double>::epsilon():exp(etai)/(opexp*opexp);
		}
	}

	Eigen::VectorXd binomial_dev_resids(const Eigen::VectorXd &y, const Eigen::VectorXd &mu, const Eigen::VectorXd &wt){
		int n = y.rows();
		int lmu = mu.rows(), lwt=wt.rows();
		Eigen::VectorXd ans = y;
		if (lmu != n && lmu != 1){
			std::string error_message = "Argument mu must be a numeric vector of length 1 or length "+std::to_string(n);
			throw std::runtime_error(error_message);
		}
		if(lwt != n && lwt != 1){
			std::string error_message = "Argument wt must be a numeric vector of length 1 or length "+std::to_string(n);
			throw std::runtime_error(error_message);
		}
		double mui, yi;
		if(lmu>1){
			for(size_t i = 0; i < n; ++i){
				mui = mu(i);
				yi = y(i);
				ans(i) = 2*wt((lwt>1)?i:0)  * (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
			}
		}
		else{
			mui = mu[0];
			for (size_t i = 0; i < n; ++i) {
			    yi = y(i);
			    ans(i) = 2 * wt((lwt > 1) ? i : 0) * (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
			}
		}
		return ans;
	}

	double binomial_dev_resids_sum(const Eigen::VectorXd &y, const Eigen::VectorXd &mu, const Eigen::VectorXd &wt){
		int n = y.rows();
		int lmu = mu.rows(), lwt=wt.rows();
		double ans = 0.0;
		if (lmu != n && lmu != 1){
			std::string error_message = "Argument mu must be a numeric vector of length 1 or length "+std::to_string(n);
			throw std::runtime_error(error_message);
		}
		if(lwt != n && lwt != 1){
			std::string error_message = "Argument wt must be a numeric vector of length 1 or length "+std::to_string(n);
			throw std::runtime_error(error_message);
		}
		double mui, yi;
		if(lmu>1){
			for(size_t i = 0; i < n; ++i){
				mui = mu(i);
				yi = y(i);
				ans+=2*wt((lwt>1)?i:0)  * (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
			}
		}
		else{
			mui = mu[0];
			for (size_t i = 0; i < n; ++i) {
				yi = y(i);
				ans+=2 * wt((lwt > 1) ? i : 0) * (y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
			}
		}
		return ans;
	}

	// This is an unsafe version of R's glm.fit
	// unsafe as in I have skipped some of the checking
	void glm(const Eigen::VectorXd &y, const Eigen::MatrixXd &x, double &p_value, double &r2, size_t max_iter, size_t thread, bool intercept){
		Eigen::setNbThreads(thread);
		/*
		Eigen::MatrixXd A;
		if(intercept){
			A=Eigen::MatrixXd::Zero(x.rows(), x.cols()+1);
			A.col(0) = Eigen::VectorXd::Constant(x.rows(),1.0);
			A.block(0,1,x.rows(), x.cols()) = x;
		}
		else A = x;
		*/
		// unfortunately, because of the algorithm where we will need to modify A, we will
		// need at least one copying
		Eigen::MatrixXd A=x;
		int nobs = y.rows();
		int nvars = A.cols();
		Eigen::VectorXd weights = Eigen::VectorXd::Ones(nobs);
		Eigen::VectorXd mustart = (weights.array()*y.array()+0.5)/(weights.array()+1);
//		Eigen::VectorXd n = Eigen::VectorXd::Ones(nobs);
//		Eigen::VectorXd m = weights.array()*y.array();
		Eigen::VectorXd eta = (mustart.array()/(1-mustart.array())).array().log();
		Eigen::VectorXd mu = logit_linkinv(eta);
		Eigen::MatrixXd A_tmp;
		double  devold=binomial_dev_resids_sum(y, mu, weights), dev=0.0;
		// Iterative reweighting
		Eigen::VectorXd varmu;
		Eigen::VectorXd mu_eta_val, z,w, good, fit, start;
		bool converge = false;
		Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr;
		qr.setThreshold(std::min(1e-7, std::numeric_limits<double>::epsilon()/1000));
		for(size_t iter = 0; iter < max_iter; ++iter){
			//varmu = (weights.array()>0).select(mu.array()*(1-mu.array()),0);
			mu_eta_val = logit_mu_eta(eta);
//			good = (weights.array()>0 && mu_eta_val.array() != 0).select(weights.array(), 0);
			good = (mu_eta_val.array() != 0).select(mu_eta_val.array(), 0);
			z = weights;
			w = weights;
			size_t i_good = 0;
//			for(size_t i_weights=0; i_weights < weights.rows(); ++i_weights){
			for(size_t i_weights=0; i_weights < good.rows(); ++i_weights){
				if(good(i_weights)>0){
					//because offset is 0, we ignore it
					z(i_good) = eta(i_weights)+(y(i_weights)-mu(i_weights))/mu_eta_val(i_weights);
					w(i_good) = std::sqrt(mu_eta_val(i_weights)*mu_eta_val(i_weights)/(mu(i_weights)*(1-mu(i_weights))));
					A.row(i_good) = A.row(i_weights);
					i_good++;
				}
			}
			z.conservativeResize(i_good);
			w.conservativeResize(i_good);
			A.conservativeResize(i_good, nvars);
			A_tmp = A;
			for(size_t i = 0; i < nvars; ++i){
				A_tmp.col(i) = A.col(i).array()*w.array();
			}
			qr.compute(A_tmp);
			start = qr.solve(Eigen::MatrixXd(z.array()*w.array()));
			if (nobs < qr.rank()){
				std::string error_message = "X matrix has rank "+std::to_string(qr.rank())+"but only "+std::to_string(nobs)+" observations";
				throw std::runtime_error(error_message);
			}
			//start = fit.transpose()*Eigen::MatrixXd(qr.colsPermutation());
			eta=A*start;
			mu = logit_linkinv(eta);
			dev = binomial_dev_resids_sum(y, mu, weights);
			// R only use 1e-8 here
			if (fabs(dev - devold)/(0.1 + fabs(dev)) < 1e-8) {
				converge =true;
				break;
			}
			else {
				devold = dev;
			}
		}
		if(!converge) throw std::runtime_error("GLM algorithm did not converge");
		Eigen::VectorXd residuals = (y.array() - mu.array())/(logit_mu_eta(eta).array());
		double rss = (residuals.array()*residuals.array()).sum();
		int sum_good = 0;
		int weight_bad = 0;
		for(size_t i = 0; i < good.rows(); ++i){
			if(good(i)!=0) sum_good++;
			if(weights(i)==0) weight_bad++;
		}
		Eigen::VectorXd wtdmu= (intercept)? Eigen::VectorXd::Constant(nobs,(weights.array()*y.array()).array().sum()/weights.array().sum()) : logit_linkinv(Eigen::VectorXd::Zero(nobs));
		double nulldev =binomial_dev_resids_sum(y, wtdmu, weights);
		int nr = std::min(sum_good, nvars);
		int n_ok = nobs - weight_bad;
		int nulldf = n_ok -intercept;
		int rank = qr.rank();
		int df_r = n_ok - rank;
		double resvar = rss/(double)df_r;
		Eigen::MatrixXd R = qr.matrixQR().topLeftCorner(rank, rank).triangularView<Eigen::Upper>();
		Eigen::VectorXd se = ((R.transpose()*R).inverse().diagonal()).array().sqrt();
		// again, we are ony interested in one of the variable
		r2 = (1.0 - std::exp((dev - nulldev)/(double)nobs))/(1 - std::exp(-nulldev/(double) nobs));
		double tvalue = start(intercept)/se(intercept);
		boost::math::normal_distribution<> dist(0,1);
		p_value = 2*boost::math::cdf(boost::math::complement(dist, fabs(tvalue)));
	}

}

