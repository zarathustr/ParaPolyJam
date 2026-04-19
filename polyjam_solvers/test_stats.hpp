#ifndef POLYJAM_SYNTHETIC_TEST_STATS_HPP_
#define POLYJAM_SYNTHETIC_TEST_STATS_HPP_

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#ifndef POLYJAM_TEST_SEEDS
#define POLYJAM_TEST_SEEDS 100
#endif

namespace polyjam_test
{

struct Stats
{
  std::vector<double> values;

  void add(double value)
  {
    values.push_back(value);
  }

  double min() const
  {
    if(values.empty()) return std::numeric_limits<double>::quiet_NaN();
    return *std::min_element(values.begin(),values.end());
  }

  double max() const
  {
    if(values.empty()) return std::numeric_limits<double>::quiet_NaN();
    return *std::max_element(values.begin(),values.end());
  }

  double mean() const
  {
    if(values.empty()) return std::numeric_limits<double>::quiet_NaN();
    return std::accumulate(values.begin(),values.end(),0.0) / static_cast<double>(values.size());
  }

  double median() const
  {
    if(values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::vector<double> sorted = values;
    std::sort(sorted.begin(),sorted.end());
    const size_t n = sorted.size();
    if(n % 2 == 1)
      return sorted[n/2];
    return 0.5 * (sorted[n/2-1] + sorted[n/2]);
  }
};

inline void printStats(const std::string & name, const Stats & stats)
{
  std::cout << std::setw(30) << std::left << name
            << " min=" << stats.min()
            << " mean=" << stats.mean()
            << " median=" << stats.median()
            << " max=" << stats.max() << std::endl;
}

inline int parseSeedCount(int argc, char ** argv)
{
  if(argc < 2)
    return POLYJAM_TEST_SEEDS;

  const int n = std::atoi(argv[1]);
  return n > 0 ? n : POLYJAM_TEST_SEEDS;
}

} // namespace polyjam_test

#endif /* POLYJAM_SYNTHETIC_TEST_STATS_HPP_ */
