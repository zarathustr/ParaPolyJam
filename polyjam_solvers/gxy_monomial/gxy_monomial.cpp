#include "gxy_monomial.hpp"

#include <cmath>

#include "gxy_monomial_core.hpp"

namespace
{

bool isFiniteVector(const Eigen::Matrix<double,6,1> & x)
{
  for(int i = 0; i < 6; ++i)
    if(!std::isfinite(x[i]))
      return false;
  return true;
}

bool isDuplicate(
    const std::vector<Eigen::Matrix<double,6,1> > & solutions,
    const Eigen::Matrix<double,6,1> & candidate)
{
  for(size_t i = 0; i < solutions.size(); ++i)
    if((solutions[i] - candidate).norm() < 1e-8)
      return true;
  return false;
}

void buildMons1(
    const Eigen::Matrix<double,6,1> & g,
    Eigen::Matrix<double,40,1> & m)
{
  const double gX1 = g[0];
  const double gX2 = g[1];
  const double gX3 = g[2];
  const double gY1 = g[3];
  const double gY2 = g[4];
  const double gY3 = g[5];

  m <<
    gX1*gY1*gY1,
    gX1*gY1*gY2,
    gX1*gY1*gY3,
    gX1*gY1,
    gX1*gY2*gY2,
    gX1*gY2*gY3,
    gX1*gY2,
    gX1*gY3*gY3,
    gX1*gY3,
    gX1,
    gX2*gY1*gY1,
    gX2*gY1*gY2,
    gX2*gY1*gY3,
    gX2*gY1,
    gX2*gY2*gY2,
    gX2*gY2*gY3,
    gX2*gY2,
    gX2*gY3*gY3,
    gX2*gY3,
    gX2,
    gX3*gY1*gY1,
    gX3*gY1*gY2,
    gX3*gY1*gY3,
    gX3*gY1,
    gX3*gY2*gY2,
    gX3*gY2*gY3,
    gX3*gY2,
    gX3*gY3*gY3,
    gX3*gY3,
    gX3,
    gY1*gY1,
    gY1*gY2,
    gY1*gY3,
    gY1,
    gY2*gY2,
    gY2*gY3,
    gY2,
    gY3*gY3,
    gY3,
    1.0;
}

void buildMons4(
    const Eigen::Matrix<double,6,1> & g,
    Eigen::Matrix<double,56,1> & m)
{
  const double gX1 = g[0];
  const double gX2 = g[1];
  const double gX3 = g[2];
  const double gY1 = g[3];
  const double gY2 = g[4];
  const double gY3 = g[5];

  m <<
    gX1*gX1*gY1,
    gX1*gX1*gY2,
    gX1*gX1*gY3,
    gX1*gX1,
    gX1*gX2*gY1,
    gX1*gX2*gY2,
    gX1*gX2*gY3,
    gX1*gX2,
    gX1*gX3*gY1,
    gX1*gX3*gY2,
    gX1*gX3*gY3,
    gX1*gX3,
    gX1*gY1,
    gX1*gY2,
    gX1*gY3,
    gX1,
    gX2*gX2*gY1,
    gX2*gX2*gY2,
    gX2*gX2*gY3,
    gX2*gX2,
    gX2*gX3*gY1,
    gX2*gX3*gY2,
    gX2*gX3*gY3,
    gX2*gX3,
    gX2*gY1,
    gX2*gY2,
    gX2*gY3,
    gX2,
    gX3*gX3*gY1,
    gX3*gX3*gY2,
    gX3*gX3*gY3,
    gX3*gX3,
    gX3*gY1,
    gX3*gY2,
    gX3*gY3,
    gX3,
    gY1*gY1*gY1,
    gY1*gY1*gY2,
    gY1*gY1*gY3,
    gY1*gY1,
    gY1*gY2*gY2,
    gY1*gY2*gY3,
    gY1*gY2,
    gY1*gY3*gY3,
    gY1*gY3,
    gY1,
    gY2*gY2*gY2,
    gY2*gY2*gY3,
    gY2*gY2,
    gY2*gY3*gY3,
    gY2*gY3,
    gY2,
    gY3*gY3*gY3,
    gY3*gY3,
    gY3,
    1.0;
}

} // namespace

void
polyjam::gxy_monomial::solve(
    const Eigen::Matrix<double,3,40> & M1,
    const Eigen::Matrix<double,3,56> & M2,
    std::vector<Eigen::Matrix<double,6,1> > & solutions)
{
  solutions.clear();

  std::vector<double> coeffs;
  coeffs.reserve(3*40 + 3*56);

  for(int r = 0; r < 3; ++r)
    for(int c = 0; c < 40; ++c)
      coeffs.push_back(M1(r,c));

  for(int r = 0; r < 3; ++r)
    for(int c = 0; c < 56; ++c)
      coeffs.push_back(M2(r,c));

  std::vector<Eigen::Matrix<double,6,1> > rawSolutions;
  polyjam::gxy_monomial_core::solve(coeffs,rawSolutions);

  for(size_t i = 0; i < rawSolutions.size(); ++i)
  {
    const Eigen::Matrix<double,6,1> candidate = rawSolutions[i];
    if(!isFiniteVector(candidate))
      continue;

    Eigen::Matrix<double,40,1> m1;
    Eigen::Matrix<double,56,1> m2;
    buildMons1(candidate,m1);
    buildMons4(candidate,m2);

    const double residual = (M1 * m1).norm() + (M2 * m2).norm();
    if(!std::isfinite(residual) || residual > 1e-5)
      continue;

    if(!isDuplicate(solutions,candidate))
      solutions.push_back(candidate);
  }
}
