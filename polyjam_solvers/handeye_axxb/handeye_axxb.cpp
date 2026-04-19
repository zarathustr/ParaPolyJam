#include "handeye_axxb.hpp"

#include <Eigen/Geometry>
#include <cmath>
#include <limits>

#include "handeye_axxb_q.hpp"

namespace
{

inline int ridx(int r, int c) { return 3*r + c; }

Eigen::Matrix3d adjugate3(const Eigen::Matrix3d & M)
{
  Eigen::Matrix3d A;

  const double a = M(0,0);
  const double b = M(0,1);
  const double c = M(0,2);
  const double d = M(1,0);
  const double e = M(1,1);
  const double f = M(1,2);
  const double g = M(2,0);
  const double h = M(2,1);
  const double i = M(2,2);

  A(0,0) = e*i - f*h;
  A(0,1) = c*h - b*i;
  A(0,2) = b*f - c*e;

  A(1,0) = f*g - d*i;
  A(1,1) = a*i - c*g;
  A(1,2) = c*d - a*f;

  A(2,0) = d*h - e*g;
  A(2,1) = b*g - a*h;
  A(2,2) = a*e - b*d;

  return A;
}

Eigen::Matrix<double,3,9> translationDesign(const Eigen::Vector3d & tB)
{
  Eigen::Matrix<double,3,9> D = Eigen::Matrix<double,3,9>::Zero();
  for(int r = 0; r < 3; ++r)
    for(int c = 0; c < 3; ++c)
      D(r,ridx(r,c)) = -tB[c];
  return D;
}

bool buildReducedObjectiveCoefficients(
    const std::vector<Eigen::Matrix4d> & A,
    const std::vector<Eigen::Matrix4d> & B,
    std::vector<double> & coeffs)
{
  coeffs.clear();
  coeffs.reserve(A.size() * 75);

  Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
  Eigen::Vector3d gAlpha = Eigen::Vector3d::Zero();
  Eigen::Matrix<double,3,9> G = Eigen::Matrix<double,3,9>::Zero();

  for(size_t i = 0; i < A.size(); ++i)
  {
    const Eigen::Matrix3d RA = A[i].block<3,3>(0,0);
    const Eigen::Vector3d tA = A[i].block<3,1>(0,3);
    const Eigen::Vector3d tB = B[i].block<3,1>(0,3);

    const Eigen::Matrix3d C = RA - Eigen::Matrix3d::Identity();
    const Eigen::Matrix<double,3,9> D = translationDesign(tB);

    H.noalias() += C.transpose() * C;
    gAlpha.noalias() += C.transpose() * tA;
    G.noalias() += C.transpose() * D;
  }

  const double detH = H.determinant();
  if(!std::isfinite(detH) || std::fabs(detH) < 1e-12)
    return false;

  const Eigen::Matrix3d adjH = adjugate3(H);
  const Eigen::Vector3d yAlpha = -adjH * gAlpha;
  const Eigen::Matrix<double,3,9> Y = -adjH * G;

  for(size_t i = 0; i < A.size(); ++i)
  {
    const Eigen::Matrix3d RA = A[i].block<3,3>(0,0);
    const Eigen::Matrix3d RB = B[i].block<3,3>(0,0);
    const Eigen::Vector3d tA = A[i].block<3,1>(0,3);
    const Eigen::Vector3d tB = B[i].block<3,1>(0,3);

    const Eigen::Matrix3d C = RA - Eigen::Matrix3d::Identity();
    const Eigen::Matrix<double,3,9> D = translationDesign(tB);

    // det(H) * vec(R_A R(q) - R(q) R_B), using exactly the same structurally
    // non-zero coefficient order as the Polyjam generator.
    for(int a = 0; a < 3; ++a)
    {
      for(int b = 0; b < 3; ++b)
      {
        for(int r = 0; r < 3; ++r)
        {
          for(int c = 0; c < 3; ++c)
          {
            if(c != b && r != a)
              continue;

            double value = 0.0;
            if(c == b)
              value += RA(a,r);
            if(r == a)
              value -= RB(c,b);

            coeffs.push_back(detH * value);
          }
        }
      }
    }

    const Eigen::Vector3d transAlpha = C * yAlpha + detH * tA;
    const Eigen::Matrix<double,3,9> transBeta = C * Y + detH * D;

    for(int k = 0; k < 3; ++k)
    {
      coeffs.push_back(transAlpha[k]);
      for(int m = 0; m < 9; ++m)
        coeffs.push_back(transBeta(k,m));
    }
  }

  for(size_t i = 0; i < coeffs.size(); ++i)
  {
    if(!std::isfinite(coeffs[i]))
      return false;
  }

  return true;
}

Eigen::Matrix3d rotationFromQuaternionSolution(
    const Eigen::Matrix<double,4,1> & q_solution,
    Eigen::Matrix<double,4,1> & q_unit,
    bool & valid)
{
  valid = false;

  const double n = q_solution.norm();
  if(!std::isfinite(n) || n < 1e-12)
    return Eigen::Matrix3d::Identity();

  q_unit = q_solution / n;

  // Keep the returned quaternion sign deterministic. q and -q are the same
  // rotation, so this only affects presentation and comparison.
  if(q_unit[0] < 0.0)
    q_unit = -q_unit;

  Eigen::Quaterniond q(q_unit[0],q_unit[1],q_unit[2],q_unit[3]);
  q.normalize();
  valid = true;
  return q.toRotationMatrix();
}

bool reconstructTranslation(
    const std::vector<Eigen::Matrix4d> & A,
    const std::vector<Eigen::Matrix4d> & B,
    const Eigen::Matrix3d & R,
    Eigen::Vector3d & t)
{
  Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
  Eigen::Vector3d g = Eigen::Vector3d::Zero();

  for(size_t i = 0; i < A.size(); ++i)
  {
    const Eigen::Matrix3d RA = A[i].block<3,3>(0,0);
    const Eigen::Vector3d tA = A[i].block<3,1>(0,3);
    const Eigen::Vector3d tB = B[i].block<3,1>(0,3);

    const Eigen::Matrix3d C = RA - Eigen::Matrix3d::Identity();
    const Eigen::Vector3d d = tA - R * tB;

    H.noalias() += C.transpose() * C;
    g.noalias() += C.transpose() * d;
  }

  Eigen::FullPivLU<Eigen::Matrix3d> lu(H);
  lu.setThreshold(1e-12);
  if(!lu.isInvertible())
    return false;

  t = -lu.solve(g);
  return t.allFinite();
}

} // namespace

void
polyjam::handeye_axxb::solve(
    std::vector<Eigen::Matrix4d> & A,
    std::vector<Eigen::Matrix4d> & B,
    std::vector<Eigen::Matrix<double,7,1> > & solutions)
{
  solutions.clear();

  if(A.size() != B.size() || A.empty())
    return;

  std::vector<double> coeffs;
  if(!buildReducedObjectiveCoefficients(A,B,coeffs))
    return;

  std::vector<Eigen::Matrix<double,4,1> > quaternion_solutions;
  polyjam::handeye_axxb_q::solve(coeffs,quaternion_solutions);

  for(size_t i = 0; i < quaternion_solutions.size(); ++i)
  {
    Eigen::Matrix<double,4,1> q_unit;
    bool valid_quaternion = false;
    const Eigen::Matrix3d R = rotationFromQuaternionSolution(
        quaternion_solutions[i],q_unit,valid_quaternion);
    if(!valid_quaternion)
      continue;

    Eigen::Vector3d t;
    if(!reconstructTranslation(A,B,R,t))
      continue;

    Eigen::Matrix<double,7,1> solution;
    solution.template block<4,1>(0,0) = q_unit;
    solution.template block<3,1>(4,0) = t;
    solutions.push_back(solution);
  }
}
