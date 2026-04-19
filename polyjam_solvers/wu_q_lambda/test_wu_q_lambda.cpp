#include <Eigen/Eigen>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "wu_q_lambda.hpp"

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

    // q0 is used when constructing a consistent Q, and qz is Polyjam's
    // default action variable.  Keep both away from zero in the test loop.
    if(std::fabs(q[0]) > 0.15 && std::fabs(q[3]) > 0.15)
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
    const double lambda,
    const Eigen::Matrix<double,4,1> & q,
    std::mt19937 & rng)
{
  const Eigen::Matrix<double,64,1> u = cubicQuaternionMonomials(q);
  Eigen::Matrix4d Q;

  for(int r = 0; r < 4; ++r)
  {
    for(int c = 1; c < 4; ++c)
      Q(r,c) = uniform(rng,1.0);

    const double target = W.row(r).dot(u) + lambda * q[r];
    double known = 0.0;
    for(int c = 1; c < 4; ++c)
      known += Q(r,c) * q[c];

    Q(r,0) = (target - known) / q[0];
  }

  return Q;
}

double systemResidual(
    const Eigen::Matrix<double,4,64> & W,
    const Eigen::Matrix4d & Q,
    const Eigen::Matrix<double,5,1> & x)
{
  const double lambda = x[0];
  Eigen::Matrix<double,4,1> q = x.template block<4,1>(1,0);
  const double n = q.norm();
  if(!std::isfinite(n) || n < 1e-12)
    return std::numeric_limits<double>::infinity();
  q /= n;

  const Eigen::Matrix<double,64,1> u = cubicQuaternionMonomials(q);
  const Eigen::Matrix<double,4,1> r = W * u - (Q - lambda * Eigen::Matrix4d::Identity()) * q;
  return r.norm() + std::fabs(q.squaredNorm() - 1.0);
}

double signInvariantQuaternionError(
    const Eigen::Matrix<double,4,1> & a,
    const Eigen::Matrix<double,4,1> & b)
{
  return std::min((a-b).norm(),(a+b).norm());
}

struct TrialResult
{
  double residual;
  double lambda_error;
  double quaternion_error;
  double time_ms;
  int solution_count;
  bool ok;
};

TrialResult runTrial(int seed)
{
  std::mt19937 rng(static_cast<unsigned int>(seed));

  Eigen::Matrix<double,4,1> q_gt = randomUnitQuaternion(rng);
  const double lambda_gt = uniform(rng,0.75);

  Eigen::Matrix<double,4,64> W;
  for(int r = 0; r < 4; ++r)
    for(int c = 0; c < 64; ++c)
      W(r,c) = uniform(rng,0.6);

  const Eigen::Matrix4d Q = buildConsistentQ(W,lambda_gt,q_gt,rng);

  std::vector<Eigen::Matrix<double,5,1> > solutions;
  const std::chrono::high_resolution_clock::time_point t0 =
      std::chrono::high_resolution_clock::now();
  polyjam::wu_q_lambda::solve(W,Q,solutions);
  const std::chrono::high_resolution_clock::time_point t1 =
      std::chrono::high_resolution_clock::now();

  const double timeMs = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(t1-t0).count();

  double bestResidual = std::numeric_limits<double>::infinity();
  double bestLambdaError = std::numeric_limits<double>::infinity();
  double bestQuaternionError = std::numeric_limits<double>::infinity();

  for(size_t i = 0; i < solutions.size(); ++i)
  {
    Eigen::Matrix<double,4,1> q = solutions[i].template block<4,1>(1,0);
    const double n = q.norm();
    if(!std::isfinite(n) || n < 1e-12)
      continue;
    q /= n;

    const double residual = systemResidual(W,Q,solutions[i]);
    const double lambdaError = std::fabs(solutions[i][0] - lambda_gt);
    const double qError = signInvariantQuaternionError(q,q_gt);

    const double score = qError + lambdaError + residual;
    const double bestScore = bestQuaternionError + bestLambdaError + bestResidual;
    if(score < bestScore)
    {
      bestResidual = residual;
      bestQuaternionError = qError;
      bestLambdaError = lambdaError;
    }
  }

  TrialResult result;
  result.residual = bestResidual;
  result.lambda_error = bestLambdaError;
  result.quaternion_error = bestQuaternionError;
  result.time_ms = timeMs;
  result.solution_count = static_cast<int>(solutions.size());
  result.ok = !solutions.empty() &&
              bestResidual <= 1e-6 &&
              bestQuaternionError <= 1e-5 &&
              bestLambdaError <= 1e-5;
  return result;
}

} // namespace

int main(int argc, char ** argv)
{
  const int numberSeeds = polyjam_test::parseSeedCount(argc,argv);

  polyjam_test::Stats residualStats;
  polyjam_test::Stats lambdaStats;
  polyjam_test::Stats quaternionStats;
  polyjam_test::Stats timeStats;
  polyjam_test::Stats solutionCountStats;

  int successCount = 0;
  for(int seed = 0; seed < numberSeeds; ++seed)
  {
    const TrialResult result = runTrial(seed);
    residualStats.add(result.residual);
    lambdaStats.add(result.lambda_error);
    quaternionStats.add(result.quaternion_error);
    timeStats.add(result.time_ms);
    solutionCountStats.add(static_cast<double>(result.solution_count));
    if(result.ok)
      ++successCount;
  }

  std::cout << "wu_q_lambda synthetic statistics over "
            << numberSeeds << " RNG seed(s)." << std::endl;
  std::cout << "Successful seeds: " << successCount << " / " << numberSeeds << std::endl;
  polyjam_test::printStats("best residual",residualStats);
  polyjam_test::printStats("lambda error",lambdaStats);
  polyjam_test::printStats("quaternion error",quaternionStats);
  polyjam_test::printStats("solve time [ms]",timeStats);
  polyjam_test::printStats("solution count",solutionCountStats);

  if(successCount != numberSeeds)
  {
    std::cerr << "wu_q_lambda statistical synthetic test failed." << std::endl;
    return 1;
  }

  std::cout << "wu_q_lambda statistical synthetic test passed." << std::endl;
  return 0;
}
