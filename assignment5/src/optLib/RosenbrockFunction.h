#pragma once

#include "ObjectiveFunction.h"

class RosenbrockFunction : public ObjectiveFunction {
public:

    RosenbrockFunction() {
		a = 1; b = 100;
    }

    virtual double computeValue(const VectorXd& x) {

		// Ex 1.1
		// return f(x)
		double xV = x[0];
		double yV = x[1];
		return (a - xV) * (a - xV) + b * (yV - (xV * xV)) * (yV - (xV * xV));
	}

    virtual void addGradientTo(VectorXd& grad, const VectorXd& x) {

		// Ex 1.1
		// write df/dx in `grad`

		double xV = x[0];
		double yV = x[1];
		grad[0] += 2 * (-a + xV + 2 * b * xV * (xV * xV - yV));
		grad[1] += 2 * b * (yV - xV * xV);

    }

	virtual void addHessianEntriesTo(std::vector<Tripletd>& hessianEntries, const VectorXd& x) {


		// Ex 1.2
		// write d^2f/dx^2 in `hessianEntries`
		double xV = x[0];
		double yV = x[1];
		hessianEntries.push_back(Tripletd(0, 0, -4 * b * (yV - xV * xV) + 8 * b * xV * xV + 2));
		hessianEntries.push_back(Tripletd(0, 1, -4 * b * xV));
		hessianEntries.push_back(Tripletd(1, 0, -4 * b * xV));
		hessianEntries.push_back(Tripletd(1, 1, 2 * b));
	}

    double a, b;
};
