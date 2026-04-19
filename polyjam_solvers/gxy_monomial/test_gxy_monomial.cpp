#include <Eigen/Eigen>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <vector>

#include "gxy_monomial.hpp"

#ifndef POLYJAM_TEST_SEEDS
#define POLYJAM_TEST_SEEDS 100
#endif

namespace
{

struct Stats
{
  std::vector<double> values;

  void add(double v)
  {
    if(std::isfinite(v))
      values.push_back(v);
  }

  double min() const
  {
    return values.empty() ? std::numeric_limits<double>::quiet_NaN()
                          : *std::min_element(values.begin(),values.end());
  }

  double max() const
  {
    return values.empty() ? std::numeric_limits<double>::quiet_NaN()
                          : *std::max_element(values.begin(),values.end());
  }

  double mean() const
  {
    if(values.empty())
      return std::numeric_limits<double>::quiet_NaN();
    double sum = 0.0;
    for(size_t i = 0; i < values.size(); ++i)
      sum += values[i];
    return sum / static_cast<double>(values.size());
  }

  double median() const
  {
    if(values.empty())
      return std::numeric_limits<double>::quiet_NaN();
    std::vector<double> sorted = values;
    std::sort(sorted.begin(),sorted.end());
    const size_t n = sorted.size();
    if(n % 2 == 1)
      return sorted[n/2];
    return 0.5 * (sorted[n/2-1] + sorted[n/2]);
  }
};

void printStats(const std::string & name, const Stats & s)
{
  std::cout << std::setw(24) << name
            << " | min=" << s.min()
            << " max=" << s.max()
            << " mean=" << s.mean()
            << " median=" << s.median()
            << std::endl;
}

Eigen::Matrix<double,6,1> randomGroundTruth(std::mt19937 & rng)
{
  std::uniform_real_distribution<double> dist(-0.8,0.8);
  for(;;)
  {
    Eigen::Matrix<double,6,1> g;
    for(int i = 0; i < 6; ++i)
      g[i] = dist(rng);

    bool good = true;
    for(int i = 0; i < 6; ++i)
      good = good && (std::fabs(g[i]) > 0.10);

    if(good)
      return g;
  }
}

double uniform(std::mt19937 & rng, double radius)
{
  std::uniform_real_distribution<double> dist(-radius,radius);
  return dist(rng);
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

void buildConsistentMeasurements(
    const Eigen::Matrix<double,6,1> & g,
    std::mt19937 & rng,
    Eigen::Matrix<double,3,40> & M1,
    Eigen::Matrix<double,3,56> & M2)
{
  Eigen::Matrix<double,40,1> m1;
  Eigen::Matrix<double,56,1> m2;
  buildMons1(g,m1);
  buildMons4(g,m2);

  for(int r = 0; r < 3; ++r)
  {
    double sum = 0.0;
    for(int c = 0; c < 39; ++c)
    {
      M1(r,c) = uniform(rng,1.0);
      sum += M1(r,c) * m1[c];
    }
    M1(r,39) = -sum;
  }

  for(int r = 0; r < 3; ++r)
  {
    double sum = 0.0;
    for(int c = 0; c < 55; ++c)
    {
      M2(r,c) = uniform(rng,1.0);
      sum += M2(r,c) * m2[c];
    }
    M2(r,55) = -sum;
  }
}

double systemResidual(
    const Eigen::Matrix<double,3,40> & M1,
    const Eigen::Matrix<double,3,56> & M2,
    const Eigen::Matrix<double,6,1> & g)
{
  Eigen::Matrix<double,40,1> m1;
  Eigen::Matrix<double,56,1> m2;
  buildMons1(g,m1);
  buildMons4(g,m2);
  return (M1 * m1).norm() + (M2 * m2).norm();
}

struct TrialResult
{
  double residual;
  double error;
  double timeMs;
  int solutionCount;
  bool success;
};

TrialResult runTrial(int seed)
{
  std::mt19937 rng(seed);

  const Eigen::Matrix<double,6,1> g_gt = randomGroundTruth(rng);

  Eigen::Matrix<double,3,40> M1;
  Eigen::Matrix<double,3,56> M2;
  buildConsistentMeasurements(g_gt,rng,M1,M2);

  std::vector<Eigen::Matrix<double,6,1> > solutions;

  const std::chrono::high_resolution_clock::time_point t0 =
      std::chrono::high_resolution_clock::now();
  polyjam::gxy_monomial::solve(M1,M2,solutions);
  const std::chrono::high_resolution_clock::time_point t1 =
      std::chrono::high_resolution_clock::now();

  double bestResidual = std::numeric_limits<double>::infinity();
  double bestError = std::numeric_limits<double>::infinity();

  for(size_t i = 0; i < solutions.size(); ++i)
  {
    const double residual = systemResidual(M1,M2,solutions[i]);
    const double error = (solutions[i] - g_gt).norm();

    if(error < bestError)
    {
      bestError = error;
      bestResidual = residual;
    }
  }

  TrialResult result;
  result.residual = bestResidual;
  result.error = bestError;
  result.timeMs = std::chrono::duration<double,std::milli>(t1-t0).count();
  result.solutionCount = static_cast<int>(solutions.size());
  result.success = bestResidual <= 1e-6 && bestError <= 1e-5;
  return result;
}

int parsePositiveInt(const char * text, int fallback)
{
  if(!text)
    return fallback;
  const int value = std::atoi(text);
  return value > 0 ? value : fallback;
}

} // namespace

int main(int argc, char ** argv)
{
  const int seedCount = argc > 1 ? parsePositiveInt(argv[1],POLYJAM_TEST_SEEDS) : POLYJAM_TEST_SEEDS;
  const int firstSeed = argc > 2 ? parsePositiveInt(argv[2],0) : 0;

  Stats residualStats;
  Stats errorStats;
  Stats timeStats;
  Stats solutionCountStats;

  int failures = 0;
  for(int i = 0; i < seedCount; ++i)
  {
    const int seed = firstSeed + i;
    const TrialResult r = runTrial(seed);

    residualStats.add(r.residual);
    errorStats.add(r.error);
    timeStats.add(r.timeMs);
    solutionCountStats.add(static_cast<double>(r.solutionCount));

    if(!r.success)
    {
      ++failures;
      std::cout << "Failed seed " << seed
                << " | residual=" << r.residual
                << " | error=" << r.error
                << " | solutions=" << r.solutionCount
                << " | time_ms=" << r.timeMs << std::endl;
    }
  }

  std::cout << "gxy_monomial synthetic statistics over " << seedCount
            << " seed(s), first_seed=" << firstSeed << std::endl;
  printStats("best_residual",residualStats);
  printStats("best_g_error",errorStats);
  printStats("solve_time_ms",timeStats);
  printStats("solution_count",solutionCountStats);
  std::cout << "failures: " << failures << std::endl;

  if(failures != 0)
  {
    std::cerr << "gxy_monomial statistical synthetic test failed." << std::endl;
    return 1;
  }

  std::cout << "gxy_monomial statistical synthetic test passed." << std::endl;
  return 0;
}
