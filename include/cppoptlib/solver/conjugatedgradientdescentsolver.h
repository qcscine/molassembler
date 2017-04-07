// CppNumericalSolver
#ifndef CONJUGATEDGRADIENTDESCENTSOLVER_H_
#define CONJUGATEDGRADIENTDESCENTSOLVER_H_

#include <Eigen/Core>
#include "isolver.h"
#include "../linesearch/armijo.h"

namespace cppoptlib {

template<typename ProblemType>
class ConjugatedGradientDescentSolver : public ISolver<ProblemType, 1> {

 public:
  using Superclass = ISolver<ProblemType, 1>;
  using typename Superclass::Scalar;
  using typename Superclass::TVector;

  /**
   * @brief minimize
   * @details [long description]
   *
   * @param objFunc [description]
   */
  void minimize(ProblemType &objFunc, TVector &x0) {
    TVector x_old(x0.size());
    Scalar value_old, value = objFunc(x0);

    TVector grad(x0.rows());
    TVector grad_old(x0.rows());
    TVector Si(x0.rows());
    TVector Si_old(x0.rows());

    this->m_current.reset();
    do {
      if(this->m_stop.xDelta > 0) x_old = x0;
      value_old = value;

      objFunc.gradient(x0, grad);

      if (this->m_current.iterations == 0) {
        Si = -grad;
      } else {
        const double beta = grad.dot(grad) / (grad_old.dot(grad_old));
        Si = -grad + beta * Si_old;
      }

      const double rate = Armijo<ProblemType, 1>::linesearch(x0, Si, objFunc);

      x0 = x0 + rate * Si;

      grad_old = grad;
      Si_old = Si;

      value = objFunc(x0);

      this->m_current.iterations += 1;
      if(this->m_stop.xDelta > 0) this->m_current.xDelta = (x0 - x_old).norm();
      this->m_current.fDelta = fabs(value - value_old);

      // TODO Why the l∞ norm? That's the maximum absolute value in the vector
      this->m_current.gradNorm = grad.template lpNorm<Eigen::Infinity>();

      this->m_status = checkConvergence(this->m_stop, this->m_current);
    } while (objFunc.callback(this->m_current, x0) && (this->m_status == Status::Continue) );

  }

  struct StepResultType {
    Scalar value, gradientNorm;
    // Supplying positions is unnecessary -> x0
    TVector negativeGradient;
    Status status;
  };

  StepResultType step(ProblemType &objFunc, TVector &x0) {
    StepResultType returnStruct;

    // Store current value, old one
    TVector x_old(x0.size());
    x_old = x0;
    Scalar value_old = objFunc(x0);

    // Make a gradient (and one to modify) 
    returnStruct.negativeGradient.resize(x0.rows());

    // Get a gradient (not negative yet)
    objFunc.gradient(x0, returnStruct.negativeGradient);

    // Negate the gradient (this is the search direction we want to follow)
    returnStruct.negativeGradient *= -1;

    // Calculate rate of position propagation
    const double rate = Armijo<ProblemType, 1>::linesearch(
      x0,
      returnStruct.negativeGradient,
      objFunc
    );

    // Propagate positions
    x0 = x0 + rate * returnStruct.negativeGradient;

    // Calculate the new value of the function
    returnStruct.value = objFunc(x0);

    // Set iterations
    this->m_current.iterations = 1;

    // Calculate convergence criteria
    if(this->m_stop.xDelta > 0) this->m_current.xDelta = (x0 - x_old).norm();
    this->m_current.fDelta = fabs(returnStruct.value - value_old);

    // TODO Why the l∞ norm? That's the maximum absolute value in the vector
    this->m_current.gradNorm = returnStruct.negativeGradient.template lpNorm<Eigen::Infinity>();
    returnStruct.gradientNorm = this -> m_current.gradNorm;

    // Check if converged
    this->m_status = checkConvergence(this->m_stop, this->m_current);
    returnStruct.status = this->m_status;

    return returnStruct;
  }

};

} /* namespace cppoptlib */

#endif /* CONJUGATEDGRADIENTDESCENTSOLVER_H_ */
