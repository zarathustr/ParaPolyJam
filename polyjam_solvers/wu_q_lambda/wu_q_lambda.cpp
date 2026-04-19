#include "wu_q_lambda.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "wu_q_lambda_core.hpp"

namespace
{

inline int cubicIndex(int a, int b, int c)
{
  return (a * 4 + b) * 4 + c;
}

Eigen::Matrix<double,64,1> cubicQuaternionMonomials(
    const Eigen::Matrix<double,4,1> & q)
{
  Eigen::Matrix<double,64,1> u;
  for(int a = 0; a < 4; ++a)
    for(int b = 0; b < 4; ++b)
      for(int c = 0; c < 4; ++c)
        u[cubicIndex(a,b,c)] = q[a] * q[b] * q[c];
  return u;
}

bool isFiniteVector(const Eigen::Matrix<double,5,1> & x)
{
  for(int i = 0; i < 5; ++i)
    if(!std::isfinite(x[i]))
      return false;
  return true;
}

double polynomialResidual(
    const Eigen::Matrix<double,4,64> & W,
    const Eigen::Matrix4d & Q,
    const double lambda,
    const Eigen::Matrix<double,4,1> & q)
{
  const Eigen::Matrix<double,64,1> u = cubicQuaternionMonomials(q);
  const Eigen::Matrix<double,4,1> r = W * u - (Q - lambda * Eigen::Matrix4d::Identity()) * q;
  return r.norm() + std::fabs(q.squaredNorm() - 1.0);
}

bool isDuplicate(
    const std::vector<Eigen::Matrix<double,5,1> > & solutions,
    const Eigen::Matrix<double,5,1> & candidate)
{
  for(size_t i = 0; i < solutions.size(); ++i)
  {
    const Eigen::Matrix<double,4,1> qi = solutions[i].template block<4,1>(1,0);
    const Eigen::Matrix<double,4,1> qc = candidate.template block<4,1>(1,0);

    const double sameSign = (qi - qc).norm();
    const double lambdaDiff = std::fabs(solutions[i][0] - candidate[0]);

    if(lambdaDiff < 1e-8 && sameSign < 1e-8)
      return true;
  }
  return false;
}

} // namespace

void
polyjam::wu_q_lambda::solve(
    const Eigen::Matrix<double,4,64> & W,
    const Eigen::Matrix4d & Q,
    std::vector<Eigen::Matrix<double,5,1> > & solutions)
{
  solutions.clear();

  std::vector<double> coeffs;
  coeffs.reserve(272);

  for(int r = 0; r < 4; ++r)
    for(int c = 0; c < 64; ++c)
      coeffs.push_back(W(r,c));

  for(int r = 0; r < 4; ++r)
    for(int c = 0; c < 4; ++c)
      coeffs.push_back(Q(r,c));

  std::vector<Eigen::Matrix<double,5,1> > rawSolutions;
  polyjam::wu_q_lambda_core::solve(coeffs,rawSolutions);

  for(size_t i = 0; i < rawSolutions.size(); ++i)
  {
    if(!isFiniteVector(rawSolutions[i]))
      continue;

    Eigen::Matrix<double,5,1> candidate = rawSolutions[i];
    Eigen::Matrix<double,4,1> q = candidate.template block<4,1>(1,0);

    const double qNorm = q.norm();
    if(!std::isfinite(qNorm) || qNorm < 1e-12)
      continue;

    q /= qNorm;
    candidate.template block<4,1>(1,0) = q;

    const double residual = polynomialResidual(W,Q,candidate[0],q);
    if(!std::isfinite(residual) || residual > 1e-5)
      continue;

    if(!isDuplicate(solutions,candidate))
      solutions.push_back(candidate);
  }
}
