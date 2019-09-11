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


#ifndef PRSICE_REGRESSION_H_
#define PRSICE_REGRESSION_H_

#include "dcdflib.h"
#include "fastlm.hpp"
#include "glm.hpp"
#include "misc.hpp"
#include <Eigen/Dense>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <stdexcept>
namespace Regression
{
void glm(const Eigen::VectorXd& y, const Eigen::MatrixXd& x, double& p_value,
         double& r2, double& coeff, double& standard_error, size_t thread = 1);
void fastLm(const Eigen::VectorXd& y, const Eigen::MatrixXd& X, double& p_value,
            double& r2, double& r2_adjust, double& coeff,
            double& standard_error, size_t thread, bool intercept,
            int type = 0);
}

#endif /* PRSICE_REGRESSION_H_ */
