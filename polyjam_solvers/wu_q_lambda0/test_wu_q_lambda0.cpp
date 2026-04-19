#include <Eigen/Eigen>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "wu_q_lambda0.hpp"

#ifndef POLYJAM_TEST_SEEDS
#define POLYJAM_TEST_SEEDS 100
#endif
#include "../test_stats.hpp"

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

Eigen::Matrix<double,4,1> randomUnitQuaternion(std::mt19937 & rng)
{
  std::uniform_real_distribution<double> dist(-1.0,1.0);

  for(;;)
  {
    Eigen::Matrix<double,4,1> q;
    q << dist(rng), dist(rng), dist(rng), dist(rng);
    const double n = q.norm();
    if(n < 1e-12)
      continue;
    q /= n;

    // q0 is used when constructing a consistent Q, and qz is commonly useful
    // as an action variable.  Keep both away from zero in the test loop.
    if(std::fabs(q[0]) > 0.20 && std::fabs(q[3]) > 0.15)
      return q;
  }
}

double uniform(std::mt19937 & rng, double radius)
{
  std::uniform_real_distribution<double> dist(-radius,radius);
  return dist(rng);
}

Eigen::Matrix4d buildConsistentQ(
    const Eigen::Matrix<double,4,64> & W,
    const Eigen::Matrix<double,4,1> & q,
    std::mt19937 & rng)
{
  const Eigen::Matrix<double,64,1> u = cubicQuaternionMonomials(q);
  Eigen::Matrix4d Q;

  for(int r = 0; r < 4; ++r)
  {
    for(int c = 1; c < 4; ++c)
      Q(r,c) = uniform(rng,1.0);

    const double target = W.row(r).dot(u);
    double known = 0.0;
    for(int c = 1; c < 4; ++c)
      known += Q(r,c) * q[c];

    Q(r,0) = (target - known) / q[0];
  }

  return Q;
}

double selectedResidual(
    const Eigen::Matrix<double,4,64> & W,
    const Eigen::Matrix4d & Q,
    const std::array<int,3> & rows,
    const Eigen::Matrix<double,4,1> & qInput)
{
  Eigen::Matrix<double,4,1> q = qInput;
  const double n = q.norm();
  if(!std::isfinite(n) || n < 1e-12)
    return std::numeric_limits<double>::infinity();
  q /= n;

  const Eigen::Matrix<double,64,1> u = cubicQuaternionMonomials(q);
  double residual = std::fabs(q.squaredNorm() - 1.0);

  for(int i = 0; i < 3; ++i)
  {
    const int r = rows[i];
    residual += std::fabs(W.row(r).dot(u) - Q.row(r).dot(q));
  }

  return residual;
}

double fullResidual(
    const Eigen::Matrix<double,4,64> & W,
    const Eigen::Matrix4d & Q,
    const Eigen::Matrix<double,4,1> & qInput)
{
  Eigen::Matrix<double,4,1> q = qInput;
  const double n = q.norm();
  if(!std::isfinite(n) || n < 1e-12)
    return std::numeric_limits<double>::infinity();
  q /= n;

  const Eigen::Matrix<double,64,1> u = cubicQuaternionMonomials(q);
  const Eigen::Matrix<double,4,1> r = W * u - Q * q;
  return r.norm() + std::fabs(q.squaredNorm() - 1.0);
}

double signInvariantQuaternionError(
    const Eigen::Matrix<double,4,1> & a,
    const Eigen::Matrix<double,4,1> & b)
{
  return std::min((a-b).norm(),(a+b).norm());
}

struct SelectionResult
{
  double selected_residual;
  double full_residual;
  double quaternion_error;
  double time_ms;
  int solution_count;
  bool ok;
};

SelectionResult runSelection(
    const Eigen::Matrix<double,4,64> & W,
    const Eigen::Matrix4d & Q,
    const Eigen::Matrix<double,4,1> & q_gt,
    const std::array<int,3> & rows)
{
  std::vector<Eigen::Matrix<double,4,1> > solutions;
  const std::chrono::high_resolution_clock::time_point t0 =
      std::chrono::high_resolution_clock::now();
  polyjam::wu_q_lambda0::solve(W,Q,rows,solutions);
  const std::chrono::high_resolution_clock::time_point t1 =
      std::chrono::high_resolution_clock::now();

  const double timeMs = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(t1-t0).count();

  double bestSelectedResidual = std::numeric_limits<double>::infinity();
  double bestFullResidual = std::numeric_limits<double>::infinity();
  double bestQuaternionError = std::numeric_limits<double>::infinity();

  for(size_t i = 0; i < solutions.size(); ++i)
  {
    Eigen::Matrix<double,4,1> q = solutions[i];
    const double n = q.norm();
    if(!std::isfinite(n) || n < 1e-12)
      continue;
    q /= n;

    const double sr = selectedResidual(W,Q,rows,q);
    const double fr = fullResidual(W,Q,q);
    const double qe = signInvariantQuaternionError(q,q_gt);

    const double score = qe + sr + fr;
    const double bestScore = bestQuaternionError + bestSelectedResidual + bestFullResidual;
    if(score < bestScore)
    {
      bestQuaternionError = qe;
      bestSelectedResidual = sr;
      bestFullResidual = fr;
    }
  }

  SelectionResult result;
  result.selected_residual = bestSelectedResidual;
  result.full_residual = bestFullResidual;
  result.quaternion_error = bestQuaternionError;
  result.time_ms = timeMs;
  result.solution_count = static_cast<int>(solutions.size());
  result.ok = !solutions.empty() &&
              bestSelectedResidual <= 1e-6 &&
              bestFullResidual <= 1e-6 &&
              bestQuaternionError <= 1e-5;
  return result;
}

} // namespace

int main(int argc, char ** argv)
{
  const int numberSeeds = polyjam_test::parseSeedCount(argc,argv);

  const std::array<int,3> rowSelections[4] = {
    std::array<int,3>{{0,1,2}},
    std::array<int,3>{{0,1,3}},
    std::array<int,3>{{0,2,3}},
    std::array<int,3>{{1,2,3}}
  };

  polyjam_test::Stats selectedResidualStats;
  polyjam_test::Stats fullResidualStats;
  polyjam_test::Stats quaternionStats;
  polyjam_test::Stats timeStats;
  polyjam_test::Stats solutionCountStats;

  int successCount = 0;
  int totalRuns = 0;

  for(int seed = 0; seed < numberSeeds; ++seed)
  {
    std::mt19937 rng(static_cast<unsigned int>(seed));

    Eigen::Matrix<double,4,1> q_gt = randomUnitQuaternion(rng);

    Eigen::Matrix<double,4,64> W;
    for(int r = 0; r < 4; ++r)
      for(int c = 0; c < 64; ++c)
        W(r,c) = uniform(rng,0.6);

    const Eigen::Matrix4d Q = buildConsistentQ(W,q_gt,rng);

    for(int i = 0; i < 4; ++i)
    {
      const SelectionResult result = runSelection(W,Q,q_gt,rowSelections[i]);
      selectedResidualStats.add(result.selected_residual);
      fullResidualStats.add(result.full_residual);
      quaternionStats.add(result.quaternion_error);
      timeStats.add(result.time_ms);
      solutionCountStats.add(static_cast<double>(result.solution_count));
      if(result.ok)
        ++successCount;
      ++totalRuns;
    }
  }

  std::cout << "wu_q_lambda0 selected-row synthetic statistics over "
            << numberSeeds << " RNG seed(s) and " << totalRuns
            << " solver call(s)." << std::endl;
  std::cout << "Successful solver calls: " << successCount << " / " << totalRuns << std::endl;
  polyjam_test::printStats("selected residual",selectedResidualStats);
  polyjam_test::printStats("full residual",fullResidualStats);
  polyjam_test::printStats("quaternion error",quaternionStats);
  polyjam_test::printStats("solve time [ms]",timeStats);
  polyjam_test::printStats("solution count",solutionCountStats);

  if(successCount != totalRuns)
  {
    std::cerr << "wu_q_lambda0 statistical synthetic test failed." << std::endl;
    return 1;
  }

  std::cout << "wu_q_lambda0 statistical synthetic test passed." << std::endl;
  return 0;
}
