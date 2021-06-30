#pragma once

#include "ObjectiveFunction.h"

#include <string>
#include <iostream>
#include <cmath>
#include <cfloat>

class GradientDescentFunctionMinimizer{
public:
	GradientDescentFunctionMinimizer(int maxIterations=100, double solveResidual=1e-5, int maxLineSearchIterations=15)
		: maxIterations(maxIterations), solveResidual(solveResidual), maxLineSearchIterations(maxLineSearchIterations){
	}

	virtual ~GradientDescentFunctionMinimizer(){}

	int getLastIterations() { return lastIterations; }

	virtual bool minimize(ObjectiveFunction *function, VectorXd &x){

		//number of parameters...
		int N = (int) x.size();
		resize(xi, N);
		resize(dx, N);
		resize(gradient, N);

		xi = x;

		bool optimizationConverged = false;

		int i=0;
		for(; i < maxIterations; i++) {
			computeSearchDirection(function, xi, dx);
			double perf = dx.norm();
			if (perf < solveResidual)
			{
				optimizationConverged = true;
				break;
			}

			doLineSearch(function, dx, xi);
		}

		lastIterations = i;

		//p now holds the parameter values at the start of the iteration...
		x = xi;

		//and done!
		return optimizationConverged;
	}

protected:
	// Since the gradient of a function gives the direction of steepest descent, all one needs to do is go in that direction...
	virtual void computeSearchDirection(ObjectiveFunction *function, const VectorXd &x, VectorXd& dx) {

		// Ex. 1.1

		dx.setZero();
		function->addGradientTo(dx, x);
		dx *= -1;
	}

	virtual void doLineSearch(ObjectiveFunction *function, const VectorXd& dx, VectorXd& xi)
	{

		// Ex. 1.1

		double alpha = 1;
		double beta = 0.5;
		bool is_inferior = false;
		VectorXd tmp_x;
		tmp_x.resize(xi.rows(), 1);
		double ref_val = function->computeValue(xi);
		double cur_alpha = alpha;

		tmp_x = xi + dx * cur_alpha;
		double tmp_val = function->computeValue(tmp_x);

		for (int i = 0; i < maxLineSearchIterations || !std::isfinite(tmp_val); i++)
		{
			if (tmp_val < ref_val && std::isfinite(tmp_val))
				break;
			cur_alpha *= beta;
			tmp_x = xi + dx * cur_alpha;
			tmp_val = function->computeValue(tmp_x);
		}

		xi = tmp_x;
	}

protected:
	double solveResidual = 1e-5;
	int maxIterations = 100;
	int maxLineSearchIterations = 15;

	VectorXd xi, dx, gradient;

	// some stats about the last time `minimize` was called
	int lastIterations = -1;
};
