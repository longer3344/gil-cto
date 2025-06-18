#pragma once

#include <Eigen/Dense>
#include <cstdint>
#include <ceres/jet.h>
#include <iomanip>
#include <vector>

namespace basalt {

// Define time types
using time_span_t = std::pair<double, double>;
using time_init_t = std::initializer_list<time_span_t>;

struct MetaData {
  virtual size_t NumParameters() const = 0;
};

template <int _N>
struct SplineSegmentMeta : public MetaData {
  static constexpr int N = _N;        ///< 样条阶数, 对于Bezier曲线，阶数和次数是一样的，但是B样条，阶数是次数加1
  static constexpr int DEG = _N - 1;  ///< 样条次数, 对于Bezier曲线，阶数和次数是一样的，但是B样条，阶数是次数加1

  double t0; ///< 样条节点的开始时间
  double dt; ///< 样条节点的时间间隔
  size_t n;  ///< 样条节点的个数

  SplineSegmentMeta(double _t0, double _dt, size_t _n = 0)
          : t0(_t0), dt(_dt), n(_n) {}

  /// @brief 样条节点的个数
  ///
  /// @return 样条节点的个数
  size_t NumParameters() const override {
    return n;
  }


  /// @brief Minimum time represented by spline
  ///
  /// @return minimum time represented by spline in seconds
  /// 实际时间要大于等于minTime
  double MinTime() const {
    return t0;
  }

  /// @brief Maximum time represented by spline
  ///
  /// @return maximum time represented by spline in nanoseconds
  /// 实际时间要小于maxTime
  double MaxTime() const {
	if (n - DEG <= 0)
	  return t0;
    return t0 + (n-DEG) * dt;
  }

  template<typename T>
  size_t PotentiallyUnsafeFloor(T x) const {
    return static_cast<size_t>(std::floor(x));
  }

  // This way of treating Jets are potentially unsafe, hence the function name
  template<typename Scalar, int N>
  size_t PotentiallyUnsafeFloor(const ceres::Jet<Scalar, N>& x) const {
    return static_cast<size_t>(ceres::floor(x.a));
  };

  /// @brief 计算当前时间的样条节点位置
  ///
  /// @param[in]  timestamp	当前时间
  /// @param[out] u         归一化时间
  /// @param[out] s         开始节点标号
  /// @return bool
  template <typename T>
  bool computeTIndex(const T& timestamp, T& u, size_t& s) const
  {  
	// 将时间调整到样条节点的有效范围内
	T t = timestamp;
    if (timestamp >= T(MaxTime()))
      t = timestamp - T(1e-9);  // 1ns
    else if(timestamp < T(MinTime()))
      t = timestamp + T(1e-9);

    if (t >= T(MinTime()) && t < T(MaxTime()))
	{
      T st = (t - T(t0)) / T(dt);     // 相对位置
      s = PotentiallyUnsafeFloor(st); // 开始节点标号
      u = st - T(s);                  // 归一化时间
      return true;
    }
	else
	{
      return false; // 不在有效时间范围内
    }

	return false;
  }
};

template <int _N>
struct SplineMeta {
  std::vector<SplineSegmentMeta<_N>> segments;

  size_t NumParameters() const  {
    size_t n = 0;
    for (auto &segment_meta : segments) {
      n += segment_meta.NumParameters();
    }
    return n;
  }

  /// @brief 计算当前时间的样条节点位置
  ///
  /// @param[in]  timestamp	当前时间
  /// @param[out] idx       开始节点位置
  /// @param[out] u			归一化时间
  /// @return bool
  template <typename T>
  bool ComputeSplineIndex(const T& timestamp, size_t& idx, T& u) const
  {
    idx = 0;
    for (auto const& seg : segments)
	{
      size_t s = 0;
      if (seg.computeTIndex(timestamp, u, s))
	  {
        idx += s;
        return true;
      } 
	  else
	  {
        idx += seg.NumParameters();
      }
    }
    std::cout << std::fixed << std::setprecision(15) << "[ComputeSplineIndex] t: " << timestamp << std::endl;
    std::cout << " not in [" << segments[0].t0 << ", " << segments[0].MaxTime() << "]" << std::endl;

    assert(timestamp >= segments[0].t0 && timestamp < segments[0].MaxTime() && "[ComputeSplineIndex] not in range");
    return false;
  }
};

}  // namespace basalt
