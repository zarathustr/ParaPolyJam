#include <Eigen/Eigen>
#include <Eigen/Geometry>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "handeye_axxb.hpp"
#include "../test_stats.hpp"

#ifndef POLYJAM_HAND_EYE_PAIRS
#define POLYJAM_HAND_EYE_PAIRS 2
#endif

#ifndef POLYJAM_TEST_SEEDS
#define POLYJAM_TEST_SEEDS 100
#endif

namespace
{

Eigen::Quaterniond randomUnitQuaternion(std::mt19937 & rng)
{
  std::uniform_real_distribution<double> dist(-1.0,1.0);

  Eigen::Vector4d qraw;
  for(;;)
  {
    qraw << dist(rng), dist(rng), dist(rng), dist(rng);
    if(qraw.norm() < 1e-12)
      continue;

    qraw.normalize();

    // The generated hand-eye solver uses q_z as its action variable. Avoid a
    // nearly zero action eigenvalue in the synthetic test loop.
    if(std::fabs(qraw[3]) >= 0.05)
      break;
  }

  if(qraw[0] < 0.0)
    qraw = -qraw;

  return Eigen::Quaterniond(qraw[0],qraw[1],qraw[2],qraw[3]);
}

Eigen::Vector3d randomVector(std::mt19937 & rng, double radius)
{
  std::uniform_real_distribution<double> dist(-radius,radius);
  return Eigen::Vector3d(dist(rng),dist(rng),dist(rng));
}

Eigen::Matrix4d makeTransform(const Eigen::Matrix3d & R, const Eigen::Vector3d & t)
{
  Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
  T.block<3,3>(0,0) = R;
  T.block<3,1>(0,3) = t;
  return T;
}

Eigen::Matrix4d inverseSE3(const Eigen::Matrix4d & T)
{
  Eigen::Matrix4d inv = Eigen::Matrix4d::Identity();
  inv.block<3,3>(0,0) = T.block<3,3>(0,0).transpose();
  inv.block<3,1>(0,3) = -inv.block<3,3>(0,0) * T.block<3,1>(0,3);
  return inv;
}

bool quaternionFromSolution(
    const Eigen::Matrix<double,7,1> & solution,
    Eigen::Quaterniond & q)
{
  q = Eigen::Quaterniond(solution[0],solution[1],solution[2],solution[3]);
  const double n = q.norm();
  if(!std::isfinite(n) || n < 1e-10)
    return false;
  q.normalize();
  return true;
}

double rotationAngle(const Eigen::Matrix3d & R)
{
  double c = 0.5 * (R.trace() - 1.0);
  if(c < -1.0) c = -1.0;
  if(c >  1.0) c =  1.0;
  return std::acos(c);
}

double axxbResidual(
    const std::vector<Eigen::Matrix4d> & A,
    const std::vector<Eigen::Matrix4d> & B,
    const Eigen::Matrix4d & X)
{
  double residual = 0.0;
  for(size_t i = 0; i < A.size(); ++i)
    residual += (A[i] * X - X * B[i]).norm();
  return residual / static_cast<double>(A.size());
}

struct TrialResult
{
  double residual;
  double rotation_error;
  double translation_error;
  double quaternion_norm_error;
  double time_ms;
  int solution_count;
  bool ok;
};

TrialResult runTrial(int seed)
{
  const int numberPairs = POLYJAM_HAND_EYE_PAIRS;
  std::mt19937 rng(static_cast<unsigned int>(seed));

  const Eigen::Quaterniond qX = randomUnitQuaternion(rng);
  const Eigen::Vector3d tX = randomVector(rng,1.0);

  const Eigen::Matrix3d RX = qX.toRotationMatrix();
  const Eigen::Matrix4d Xgt = makeTransform(RX,tX);
  const Eigen::Matrix4d XgtInv = inverseSE3(Xgt);

  std::vector<Eigen::Matrix4d> A;
  std::vector<Eigen::Matrix4d> B;
  A.reserve(numberPairs);
  B.reserve(numberPairs);

  for(int i = 0; i < numberPairs; ++i)
  {
    const Eigen::Quaterniond qB = randomUnitQuaternion(rng);
    const Eigen::Vector3d tB = randomVector(rng,1.0);
    const Eigen::Matrix4d Bi = makeTransform(qB.toRotationMatrix(),tB);
    const Eigen::Matrix4d Ai = Xgt * Bi * XgtInv;

    A.push_back(Ai);
    B.push_back(Bi);
  }

  std::vector<Eigen::Matrix<double,7,1> > solutions;
  const std::chrono::high_resolution_clock::time_point t0 =
      std::chrono::high_resolution_clock::now();
  polyjam::handeye_axxb::solve(A,B,solutions);
  const std::chrono::high_resolution_clock::time_point t1 =
      std::chrono::high_resolution_clock::now();

  const double timeMs = std::chrono::duration_cast<std::chrono::duration<double,std::milli> >(t1-t0).count();

  double bestResidual = std::numeric_limits<double>::infinity();
  double bestRotationError = std::numeric_limits<double>::infinity();
  double bestTranslationError = std::numeric_limits<double>::infinity();
  double bestQuaternionNormError = std::numeric_limits<double>::infinity();

  for(size_t i = 0; i < solutions.size(); ++i)
  {
    Eigen::Quaterniond q;
    if(!quaternionFromSolution(solutions[i],q))
      continue;

    const Eigen::Vector3d t = solutions[i].block<3,1>(4,0);
    const Eigen::Matrix3d R = q.toRotationMatrix();
    const Eigen::Matrix4d X = makeTransform(R,t);

    const double residual = axxbResidual(A,B,X);
    const double rotationError = rotationAngle(RX.transpose() * R);
    const double translationError = (t - tX).norm();
    const double quaternionNormError = std::fabs(solutions[i].block<4,1>(0,0).norm() - 1.0);

    const double score = rotationError + translationError + residual;
    const double bestScore = bestRotationError + bestTranslationError + bestResidual;
    if(score < bestScore)
    {
      bestResidual = residual;
      bestRotationError = rotationError;
      bestTranslationError = translationError;
      bestQuaternionNormError = quaternionNormError;
    }
  }

  TrialResult result;
  result.residual = bestResidual;
  result.rotation_error = bestRotationError;
  result.translation_error = bestTranslationError;
  result.quaternion_norm_error = bestQuaternionNormError;
  result.time_ms = timeMs;
  result.solution_count = static_cast<int>(solutions.size());
  result.ok = !solutions.empty() &&
              bestResidual < 1e-6 &&
              bestRotationError < 1e-5 &&
              bestTranslationError < 1e-5 &&
              bestQuaternionNormError < 1e-5;
  return result;
}

} // namespace

int main(int argc, char ** argv)
{
  const int numberSeeds = polyjam_test::parseSeedCount(argc,argv);

  polyjam_test::Stats residualStats;
  polyjam_test::Stats rotationStats;
  polyjam_test::Stats translationStats;
  polyjam_test::Stats qnormStats;
  polyjam_test::Stats timeStats;
  polyjam_test::Stats solutionCountStats;

  int successCount = 0;
  for(int seed = 0; seed < numberSeeds; ++seed)
  {
    const TrialResult result = runTrial(seed);
    residualStats.add(result.residual);
    rotationStats.add(result.rotation_error);
    translationStats.add(result.translation_error);
    qnormStats.add(result.quaternion_norm_error);
    timeStats.add(result.time_ms);
    solutionCountStats.add(static_cast<double>(result.solution_count));
    if(result.ok)
      ++successCount;
  }

  std::cout << "Synthetic quaternion hand-eye AX=XB statistics over "
            << numberSeeds << " RNG seed(s)." << std::endl;
  std::cout << "Successful seeds: " << successCount << " / " << numberSeeds << std::endl;
  polyjam_test::printStats("best AX=XB residual",residualStats);
  polyjam_test::printStats("rotation error [rad]",rotationStats);
  polyjam_test::printStats("translation error",translationStats);
  polyjam_test::printStats("quaternion norm error",qnormStats);
  polyjam_test::printStats("solve time [ms]",timeStats);
  polyjam_test::printStats("solution count",solutionCountStats);

  if(successCount != numberSeeds)
  {
    std::cerr << "Synthetic quaternion hand-eye AX=XB statistical test failed." << std::endl;
    return 1;
  }

  std::cout << "Synthetic quaternion hand-eye AX=XB statistical test passed." << std::endl;
  return 0;
}
