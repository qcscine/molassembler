#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <numeric>

namespace numeric {

unsigned average(
  const std::vector<unsigned>& values
) {
  if(values.size() == 0) return 0;
  else return (double) std::floor(
    (double) std::accumulate(
      values.begin(),
      values.end(),
      0
    ) / (double) values.size()
  );
}

unsigned population_stddev(
  const std::vector<unsigned>& values
) {
  long double average = std::accumulate(
    values.begin(),
    values.end(),
    0l
  ) / (long double) values.size();

  long double intermediate = std::accumulate(
    values.begin(),
    values.end(),
    0l,
    [&average](const long double& carry, const unsigned& value) {
      return (
        carry 
        + std::pow(
          (long double) value - average,
          2l
        )
      );
    }
  );
  return std::round(
    std::sqrt( intermediate / (long double) values.size() )
  );
}


std::vector<unsigned> map_average(
  const std::vector<
    std::vector<unsigned>
  >& list_values
) {
  std::vector<unsigned> mapped;
  for(const auto& values: list_values) {
    mapped.push_back(
      average(values)
    );
  }
  return mapped;
}

std::vector<unsigned> map_population_stddev(
  const std::vector<
    std::vector<unsigned>
  >& list_values
) {
  std::vector<unsigned> mapped;
  for(unsigned i = 0; i < list_values.size(); i++) {
    mapped.push_back(
      population_stddev(list_values[i])
    );
  }

  return mapped;
}

} // eo namespace
