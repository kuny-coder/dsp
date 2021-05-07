#ifndef DSP_DSP_H
#define DSP_DSP_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <numeric>
#include <type_traits>
#include <vector>

template <typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 42>
using Point = std::vector<T>;

template <typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 42>
T discrepancy(const std::vector<Point<T>>& points, const std::vector<std::size_t>& inputIdxes, const std::vector<std::pair<T, T>>& boundaries, int subsample) {
  const size_t D = boundaries.size();
  long double part1 = 1.0 / std::powl(3.0, D);

  std::vector<std::size_t> idxes = inputIdxes;
  if (subsample > 0 && subsample < inputIdxes.size()) {
    std::random_shuffle(idxes.begin(), idxes.end());
    idxes.resize(subsample);
  }

  long double part2 = 0.0;
  for (size_t i = 0; i < idxes.size(); ++i) {
    long double prod = 1.0;
    for (size_t k = 0; k < D; ++k) {
      const auto scaled = (points[idxes[i]][k] - boundaries[k].first) / (boundaries[k].second - boundaries[k].first);
      prod *= (1.0 - std::powl(scaled, 2));
    }
    part2 += prod;
  }
  part2 /= (std::powl(2.0, D - 1) * idxes.size());

  long double part3 = 0.0;
  for (size_t i = 0; i < idxes.size(); ++i) {
    for (size_t j = 0; j < idxes.size(); ++j) {
      long double prod = 1.0;
      for (size_t k = 0; k < D; ++k) {
        const auto scaled1 = (points[idxes[i]][k] - boundaries[k].first) / (boundaries[k].second - boundaries[k].first);
        const auto scaled2 = (points[idxes[j]][k] - boundaries[k].first) / (boundaries[k].second - boundaries[k].first);
        prod *= std::min(1.0 - scaled1, 1.0 - scaled2);
      }
      part3 += prod;
    }
  }
  part3 /= (idxes.size() * idxes.size());
  return std::sqrtl(std::max<long double>(part1 - part2 + part3, 0.0));
}

template <typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 42>
std::pair<std::pair<size_t, size_t>, std::pair<std::vector<size_t>, std::vector<size_t>>>
selectMaxDim(const std::vector<Point<T>>& points, const std::vector<std::pair<T, T>>& boundaries, const std::vector<std::size_t>& idxes, size_t splits) {
  const std::size_t D = boundaries.size();
  long double maxGap = -1.0;
  std::size_t dim, splitIdx;
  std::vector<size_t> pointIdxesL, pointIdxesR;

  std::vector<size_t> idxesL, idxesR;
  idxesL.reserve(idxes.size());
  idxesR.reserve(idxes.size());

  dim = 0;
  for (std::size_t i = 1; i < D; ++i) {
    if (boundaries[dim].second - boundaries[dim].first < boundaries[i].second - boundaries[i].first) {
      dim = i;
    }
  }

  for (std::size_t j = 1; j < splits; ++j) {
    idxesL.clear();
    idxesR.clear();
    for (std::size_t l = 0; l < idxes.size(); ++l) {
      if (points[idxes[l]][dim] < boundaries[dim].first + (boundaries[dim].second - boundaries[dim].first) * j / splits) {
        idxesL.push_back(idxes[l]);
      } else {
        idxesR.push_back(idxes[l]);
      }
    }
    long double gap = abs(j * 1.0 / splits - idxesL.size() * 1.0 / idxes.size());
    if (maxGap < gap) {
      maxGap = gap;
      splitIdx = j;
      pointIdxesL = idxesL;
      pointIdxesR = idxesR;
    }
  }
  return {{dim, splitIdx}, {pointIdxesL, pointIdxesR}};
}

template <typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 42>
std::pair<std::pair<size_t, size_t>, std::pair<std::vector<size_t>, std::vector<size_t>>>
selectMaxGap(const std::vector<Point<T>>& points, const std::vector<std::pair<T, T>>& boundaries, const std::vector<std::size_t>& idxes, size_t splits, double maxDimRatio) {
  const std::size_t D = boundaries.size();
  long double maxGap = -1.0;
  std::size_t dim, splitIdx;
  std::vector<size_t> pointIdxesL, pointIdxesR;

  std::vector<size_t> idxesL, idxesR;
  idxesL.reserve(idxes.size());
  idxesR.reserve(idxes.size());

  std::size_t minDim = 0;
  std::size_t maxDim = 0;
  for (std::size_t i = 1; i < D; ++i) {
    if (boundaries[i].second - boundaries[i].first < boundaries[minDim].second - boundaries[minDim].first) {
      minDim = i;
    }
    if (boundaries[i].second - boundaries[i].first > boundaries[maxDim].second - boundaries[maxDim].first) {
      maxDim = i;
    }
  }

  if ((boundaries[maxDim].second - boundaries[maxDim].first) / (boundaries[minDim].second - boundaries[minDim].first) > maxDimRatio) {
    return selectMaxDim(points, boundaries, idxes, splits);
  }

  for (std::size_t i = 0; i < D; ++i) {
    for (std::size_t j = 1; j < splits; ++j) {
      idxesL.clear();
      idxesR.clear();
      for (std::size_t l = 0; l < idxes.size(); ++l) {
        if (points[idxes[l]][i] < boundaries[i].first + (boundaries[i].second - boundaries[i].first) * j / splits) {
          idxesL.push_back(idxes[l]);
        } else {
          idxesR.push_back(idxes[l]);
        }
      }
      long double gap = abs(j * 1.0 / splits - idxesL.size() * 1.0 / idxes.size());
      if (maxGap < gap) {
        maxGap = gap;
        dim = i;
        splitIdx = j;
        pointIdxesL = idxesL;
        pointIdxesR = idxesR;
      }
    }
  }
  return {{dim, splitIdx}, {pointIdxesL, pointIdxesR}};
}

template <typename T, typename std::enable_if<std::is_floating_point<T>::value, int>::type = 42>
std::vector<std::pair<std::vector<std::pair<T, T>>, T>>
DSP(const std::vector<Point<T>>& points, const std::vector<std::pair<T, T>>& boundaries, size_t splits, double maxDimRatio, T theta, int subsample, std::size_t minsample) {
  const long double totalRoot = std::sqrtl(points.size());
  std::vector<std::pair<std::vector<std::pair<T, T>>, T>> ret{{boundaries, 1.0}};
  std::vector<std::vector<std::size_t>> pointsIn{std::vector<std::size_t>(points.size())};
  std::iota(pointsIn[0].begin(), pointsIn[0].end(), 0);
  while (true) {
    std::vector<std::pair<std::vector<std::pair<T, T>>, T>> subs{};
    std::vector<std::vector<std::size_t>> pointsInSubs{};
    for (size_t i = 0; i < ret.size(); ++i) {
      T ratio{};
      if (pointsIn[i].size() > minsample) {
        ratio = theta * totalRoot / pointsIn[i].size();
        if (subsample > 0) {
          ratio *= std::sqrtl(pointsIn[i].size()) / std::sqrtl(std::min<size_t>(pointsIn[i].size(), subsample));
        }
      }
      if (pointsIn[i].size() > minsample && discrepancy(points, pointsIn[i], ret[i].first, subsample) > theta * ratio) {
        const auto gap = selectMaxGap(points, ret[i].first, pointsIn[i], splits, maxDimRatio);
        const std::pair<T, T> dim = ret[i].first[gap.first.first];
        const T mass = ret[i].second;
        std::vector<std::pair<T, T>> l = ret[i].first;
        l[gap.first.first].second = dim.first + (dim.second - dim.first) / splits * gap.first.second;
        std::vector<std::pair<T, T>> r = ret[i].first;
        r[gap.first.first].first = dim.first + (dim.second - dim.first) / splits * gap.first.second;
        subs.push_back({l, gap.second.first.size() * mass / pointsIn[i].size()});
        subs.push_back({r, gap.second.second.size() * mass / pointsIn[i].size()});
        pointsInSubs.push_back(gap.second.first);
        pointsInSubs.push_back(gap.second.second);
      } else {
        subs.push_back(std::move(ret[i]));
        pointsInSubs.push_back(std::move(pointsIn[i]));
      }
    }
    if (ret.size() != subs.size()) {
      ret = subs;
      pointsIn = pointsInSubs;
    } else {
      ret = subs;
      pointsIn = pointsInSubs;
      break;
    }
  }
  return ret;
}

#endif //DSP_DSP_H
