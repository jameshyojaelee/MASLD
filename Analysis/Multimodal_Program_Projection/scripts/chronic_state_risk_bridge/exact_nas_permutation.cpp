#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static std::vector<double> ranks(const std::vector<double>& x) {
  std::vector<size_t> order(x.size());
  for (size_t i = 0; i < x.size(); ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return x[a] < x[b]; });
  std::vector<double> output(x.size());
  for (size_t left = 0; left < order.size();) {
    size_t right = left + 1;
    while (right < order.size() && x[order[right]] == x[order[left]]) ++right;
    const double rank = (static_cast<double>(left + 1) + static_cast<double>(right)) / 2.0;
    for (size_t k = left; k < right; ++k) output[order[k]] = rank;
    left = right;
  }
  return output;
}

static double correlation(const std::vector<double>& x, const std::vector<double>& y) {
  double mx = 0, my = 0;
  for (size_t i = 0; i < x.size(); ++i) { mx += x[i]; my += y[i]; }
  mx /= x.size(); my /= y.size();
  double numerator = 0, dx = 0, dy = 0;
  for (size_t i = 0; i < x.size(); ++i) {
    numerator += (x[i] - mx) * (y[i] - my);
    dx += (x[i] - mx) * (x[i] - mx); dy += (y[i] - my) * (y[i] - my);
  }
  return numerator / std::sqrt(dx * dy);
}

int main(int argc, char** argv) {
  if (argc != 3) { std::cerr << "usage: exact_nas_permutation input.tsv output.tsv\n"; return 2; }
  std::ifstream in(argv[1]); std::string line; std::getline(in, line);
  std::vector<double> state, nas;
  while (std::getline(in, line)) {
    std::stringstream stream(line); std::string left, right;
    std::getline(stream, left, '\t'); std::getline(stream, right, '\t');
    if (!left.empty() && !right.empty()) { state.push_back(std::stod(left)); nas.push_back(std::stod(right)); }
  }
  if (state.size() != 13) { std::cerr << "expected 13 primary pairs, found " << state.size() << "\n"; return 3; }
  const std::vector<double> state_rank = ranks(state);
  std::vector<double> permuted_rank = ranks(nas);
  std::sort(permuted_rank.begin(), permuted_rank.end());
  const double observed = correlation(state_rank, ranks(nas));
  unsigned long long total = 0, extreme = 0;
  do {
    const double current = correlation(state_rank, permuted_rank);
    ++total; if (std::abs(current) >= std::abs(observed) - 1e-14) ++extreme;
  } while (std::next_permutation(permuted_rank.begin(), permuted_rank.end()));
  std::ofstream out(argv[2]);
  out << "estimate\tp\tn_permutations\tn_extreme\n" << std::setprecision(17)
      << observed << '\t' << static_cast<double>(extreme) / static_cast<double>(total) << '\t'
      << total << '\t' << extreme << '\n';
  return total == 10810800ULL ? 0 : 4;
}
