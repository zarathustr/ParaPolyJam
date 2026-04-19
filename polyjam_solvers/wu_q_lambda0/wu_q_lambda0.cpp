#include "wu_q_lambda0.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "wu_q_lambda0_core.hpp"

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

bool validSelectedRows(const std::array<int,3> & rows)
{
  bool used[4] = {false,false,false,false};
  for(int i = 0; i < 3; ++i)
  {
    const int r = rows[i];
    if(r < 0 || r >= 4 || used[r])
      return false;
    used[r] = true;
  }
  return true;
}

bool isFiniteVector(const Eigen::Matrix<double,4,1> & x)
{
  for(int i = 0; i < 4; ++i)
    if(!std::isfinite(x[i]))
      return false;
  return true;
}

double selectedPolynomialResidual(
    const Eigen::Matrix<double,4,64> & W,
    const Eigen::Matrix4d & Q,
    const std::array<int,3> & selected_rows,
    const Eigen::Matrix<double,4,1> & q)
{
  const Eigen::Matrix<double,64,1> u = cubicQuaternionMonomials(q);
  double residual = std::fabs(q.squaredNorm() - 1.0);

  for(int i = 0; i < 3; ++i)
  {
    const int r = selected_rows[i];
    residual += std::fabs(W.row(r).dot(u) - Q.row(r).dot(q));
  }

  return residual;
}

bool isDuplicate(
    const std::vector<Eigen::Matrix<double,4,1> > & solutions,
    const Eigen::Matrix<double,4,1> & candidate)
{
  for(size_t i = 0; i < solutions.size(); ++i)
  {
    const double sameSign = (solutions[i] - candidate).norm();
    const double oppositeSign = (solutions[i] + candidate).norm();
    if(std::min(sameSign,oppositeSign) < 1e-8)
      return true;
  }
  return false;
}

} // namespace

void
polyjam::wu_q_lambda0::solve(
    const Eigen::Matrix<double,4,64> & W,
    const Eigen::Matrix4d & Q,
    std::vector<Eigen::Matrix<double,4,1> > & solutions)
{
  const std::array<int,3> defaultRows = {{0,1,2}};
  solve(W,Q,defaultRows,solutions);
}

void
polyjam::wu_q_lambda0::solve(
    const Eigen::Matrix<double,4,64> & W,
    const Eigen::Matrix4d & Q,
    const std::array<int,3> & selected_rows,
    std::vector<Eigen::Matrix<double,4,1> > & solutions)
{
  solutions.clear();

  if(!validSelectedRows(selected_rows))
    return;

  std::vector<double> coeffs;
  coeffs.reserve(3*64 + 3*4);

  for(int i = 0; i < 3; ++i)
  {
    const int r = selected_rows[i];
    for(int c = 0; c < 64; ++c)
      coeffs.push_back(W(r,c));
  }

  for(int i = 0; i < 3; ++i)
  {
    const int r = selected_rows[i];
    for(int c = 0; c < 4; ++c)
      coeffs.push_back(Q(r,c));
  }

  std::vector<Eigen::Matrix<double,4,1> > rawSolutions;
  polyjam::wu_q_lambda0_core::solve(coeffs,rawSolutions);

  for(size_t i = 0; i < rawSolutions.size(); ++i)
  {
    if(!isFiniteVector(rawSolutions[i]))
      continue;

    Eigen::Matrix<double,4,1> q = rawSolutions[i];

    const double qNorm = q.norm();
    if(!std::isfinite(qNorm) || qNorm < 1e-12)
      continue;

    q /= qNorm;

    // q and -q represent the same quaternion rotation and are both roots of
    // this odd polynomial system.  Return a deterministic sign.
    if(q[0] < 0.0)
      q = -q;

    const double residual = selectedPolynomialResidual(W,Q,selected_rows,q);
    if(!std::isfinite(residual) || residual > 1e-5)
      continue;

    if(!isDuplicate(solutions,q))
      solutions.push_back(q);
  }
}
