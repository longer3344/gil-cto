/*
 * OA-LICalib:
 * Observability-Aware Intrinsic and Extrinsic Calibration of LiDAR-IMU Systems
 *
 * Copyright (C) 2022 Jiajun Lv
 * Copyright (C) 2022 Xingxing Zuo
 * Copyright (C) 2022 Kewei Hu
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef GNSS_FACTOR_H
#define GNSS_FACTOR_H

#include <basalt/spline/ceres_spline_helper.h>
#include <basalt/spline/spline_segment.h>
#include <ceres/ceres.h>
#include <sensor_data/gnss_data.h>
#include <sophus/so3.hpp>

static Eigen::Matrix3d MatrixSkewSymmetric(const Eigen::Vector3d v)
{
	Eigen::Matrix3d M;
	M(0,0) =  0.0; M(0,1) =-v(2); M(0,2) = v(1);
	M(1,0) = v[2]; M(1,1) = 0.0;  M(1,2) =-v(0);
	M(2,0) =-v[1]; M(2,1) = v[0]; M(2,2) =  0.0;
	return M;
}

namespace liso {
using namespace basalt;

// 
class GNSSPositionFactor : public CeresSplineHelper<SplineOrder> {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  GNSSPositionFactor(const double timestamp, const Eigen::Vector3d &position, const Eigen::Vector3d &blarm,
                                  const SplineMeta<SplineOrder>& spline_meta, double weight)
      : timestamp_(timestamp), position_(position), blarm_(blarm), spline_meta_(spline_meta) {
	weight_v3d_(0) = weight_v3d_(1) = weight_v3d_(2) = weight;
    inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
  }

  GNSSPositionFactor(const double timestamp, const Eigen::Vector3d &position, const Eigen::Vector3d &blarm,
	  const SplineMeta<SplineOrder>& spline_meta, Eigen::Vector3d weight_v3d)
	  : timestamp_(timestamp), position_(position), blarm_(blarm), spline_meta_(spline_meta), weight_v3d_(weight_v3d) {
	  inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
  }

  template <class T>
  bool operator()(T const* const* sKnots, T* sResiduals) const {
	using Vec3T = Eigen::Matrix<T, 3, 1>;
	using Vec6T = Eigen::Matrix<T, 6, 1>;
	using SO3T = Sophus::SO3<T>;
	using Tangent = typename Sophus::SO3<T>::Tangent;

	// sKnots R_ItoW/p_IinW
	size_t R_offset; // 姿态
	size_t P_offset; // 位置
	double u;
	spline_meta_.ComputeSplineIndex(fmod(timestamp_, 604800.0), R_offset, u);
	P_offset = R_offset + spline_meta_.NumParameters();

	SO3T R_IkToG;
	CeresSplineHelper<SplineOrder>::template evaluate_lie<T, Sophus::SO3>(sKnots + R_offset, u, inv_dt_, &R_IkToG);

	Vec3T p_IkinG;
	CeresSplineHelper<SplineOrder>::template evaluate<T, 3, 0>(sKnots + P_offset, u, inv_dt_, &p_IkinG);

	Eigen::Map<Vec3T> residuals(sResiduals);
	residuals =  p_IkinG + R_IkToG * blarm_ - position_;
	//std::cout << timestamp_ << std::endl;
	//std::cout << residuals << std::endl;
	//std::cout << weight_v3d_ << std::endl;

	residuals[0] = T(1.0 / weight_v3d_(0)) * residuals[0];
	residuals[1] = T(1.0 / weight_v3d_(1)) * residuals[1];
	residuals[2] = T(1.0 / weight_v3d_(2)) * residuals[2];
	//std::cout << timestamp_ << std::endl;
	//std::cout << residuals << std::endl;

	return true;
  }

 private:
  double timestamp_;
  Eigen::Vector3d  position_;
  Eigen::Vector3d  blarm_;
  SplineMeta<SplineOrder> spline_meta_;
  Eigen::Vector3d weight_v3d_;
  double inv_dt_;
};

class GNSSVelocityFactor : public CeresSplineHelper<SplineOrder> {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	GNSSVelocityFactor(const double timestamp, const Eigen::Vector3d &velocity, const Eigen::Vector3d &blarm,
		const SplineMeta<SplineOrder>& spline_meta, double weight)
		: timestamp_(timestamp), velocity_(velocity), blarm_(blarm), spline_meta_(spline_meta) {
		weight_v3d_(0) = weight_v3d_(1) = weight_v3d_(2) = weight;
		inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
	}

	GNSSVelocityFactor(const double timestamp, const Eigen::Vector3d &velocity, const Eigen::Vector3d &blarm,
		const SplineMeta<SplineOrder>& spline_meta, Eigen::Vector3d weight_v3d)
		: timestamp_(timestamp), velocity_(velocity), blarm_(blarm), spline_meta_(spline_meta), weight_v3d_(weight_v3d) {
		
		inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
	}

	template <class T>
	bool operator()(T const* const* sKnots, T* sResiduals) const {
		using Vec3T = Eigen::Matrix<T, 3, 1>;
		using Vec6T = Eigen::Matrix<T, 6, 1>;
		using Mat3T = Eigen::Matrix<T, 3, 3>;
		using SO3T = Sophus::SO3<T>;
		using Tangent = typename Sophus::SO3<T>::Tangent;

		// sKnots R_ItoW/p_IinW
		size_t R_offset; // 姿态节点
		size_t P_offset; // 位置节点
		double u;
		spline_meta_.ComputeSplineIndex(fmod(timestamp_, 604800.0), R_offset, u);
		P_offset = R_offset + spline_meta_.NumParameters();

		// 姿态
		SO3T	R_IkToG; // Rbe
		Tangent	R_IkToG_dot; // wibb
		Vec3T	wibb;
		Mat3T	wibbx;
		CeresSplineHelper<SplineOrder>::template evaluate_lie<T, Sophus::SO3>(sKnots + R_offset, u, inv_dt_, &R_IkToG, &R_IkToG_dot);
		wibb << R_IkToG_dot(0), R_IkToG_dot(1), R_IkToG_dot(2);
		wibbx << T(0), -R_IkToG_dot(2), R_IkToG_dot(1), R_IkToG_dot(2), T(0), -R_IkToG_dot(0), -R_IkToG_dot(1), R_IkToG_dot(0), T(0);

		// 速度
		Vec3T v_IkinG;
		CeresSplineHelper<SplineOrder>::template evaluate<T, 3, 1>(sKnots + P_offset, u, inv_dt_, &v_IkinG);

		Mat3T wieex  = MatrixSkewSymmetric(Eigen::Vector3d(0.0, 0.0, 7.2921151467E-5)).template cast<T>(); // 地球自转
		Mat3T blarmx = MatrixSkewSymmetric(blarm_).template cast<T>(); // 地球自转

		Eigen::Map<Vec3T> residuals(sResiduals);
		residuals = v_IkinG - R_IkToG.matrix() * blarmx * wibb - wieex * R_IkToG.matrix() * blarm_ - velocity_;
		residuals[0] = T(1.0 / weight_v3d_(0)) * residuals[0];
		residuals[1] = T(1.0 / weight_v3d_(1)) * residuals[1];
		residuals[2] = T(1.0 / weight_v3d_(2)) * residuals[2];

		return true;
	}

private:
	double timestamp_;
	Eigen::Vector3d  velocity_;
	Eigen::Vector3d  weight_v3d_;
	Eigen::Vector3d  blarm_;
	SplineMeta<SplineOrder> spline_meta_;
	double inv_dt_;
};

}  // namespace liso

#endif
