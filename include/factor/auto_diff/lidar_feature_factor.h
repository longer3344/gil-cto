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

#ifndef AUTO_DIFF_LIDAR_FEATURE_FACTOR
#define AUTO_DIFF_LIDAR_FEATURE_FACTOR

#include <basalt/spline/ceres_spline_helper_jet.h>
#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <sensor_data/imu_data.h>
#include <sensor_data/lidar_feature.h>
#include <Eigen/Core>
#include <memory>
#include <sophus/so3.hpp>

#include "LiDARFFStream.h"

namespace liso {

using namespace basalt;

class PointFeatureFactor {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  PointFeatureFactor(const PointCorrespondence& pc,
                     const SplineMeta<SplineOrder>& spline_meta, double weight)
      : measurement_(pc), spline_meta_(spline_meta), weight_(weight) {
    inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
  }

  template <class T>
  bool operator()(T const* const* sKnots, T* sResiduals) const {
    using SO3T = Sophus::SO3<T>;
    using Vec3T = Eigen::Matrix<T, 3, 1>;

	// 时间
    T t[2];
    t[0] = T(measurement_.t_map);
    t[1] = T(measurement_.t_point);

	// 计算当前时刻的IMU位姿
    T u[2];
    size_t R_offset[2];
    size_t P_offset[2];
    spline_meta_.ComputeSplineIndex(t[0], R_offset[0], u[0]); // 参考时刻, 建图时刻
    P_offset[0] = R_offset[0] + spline_meta_.NumParameters();
 
    spline_meta_.ComputeSplineIndex(t[1], R_offset[1], u[1]); // 当前时刻
    P_offset[1] = R_offset[1] + spline_meta_.NumParameters();

    SO3T R_IkToG[2];
    CeresSplineHelperJet<T, SplineOrder>::template evaluate_lie<Sophus::SO3>(
        sKnots + R_offset[0], u[0], inv_dt_, &R_IkToG[0]); // 建图时刻的姿态
    CeresSplineHelperJet<T, SplineOrder>::template evaluate_lie<Sophus::SO3>(
        sKnots + R_offset[1], u[1], inv_dt_, &R_IkToG[1]); // 当前时刻的姿态

    Vec3T p_IkinG[2];
    CeresSplineHelperJet<T, SplineOrder>::template evaluate<3, 0>(
        sKnots + P_offset[0], u[0], inv_dt_, &p_IkinG[0]); // 建图时刻的姿态
    CeresSplineHelperJet<T, SplineOrder>::template evaluate<3, 0>(
        sKnots + P_offset[1], u[1], inv_dt_, &p_IkinG[1]); // 当前时刻的位置

	// 外参
    size_t Kont_offset = 2 * spline_meta_.NumParameters();
    Eigen::Map<SO3T const> const R_LtoI(sKnots[Kont_offset]);
    Eigen::Map<Vec3T const> const p_LinI(sKnots[Kont_offset + 1]);

	// 计算相对位姿
    SO3T R_IkToI0 = R_IkToG[0].inverse() * R_IkToG[1];
    Vec3T p_IkinI0 = R_IkToG[0].inverse() * (p_IkinG[1] - p_IkinG[0]);

    // {M} frame is coinside with {L0) frame.
    Vec3T p_Lk = measurement_.point.template cast<T>(); ///< 当前时刻雷达系下的原始点
    Vec3T p_Ik = R_LtoI * p_Lk + p_LinI;				///< 当前时刻惯导系下的原始点
    Vec3T p_I0 = R_IkToI0 * p_Ik + p_IkinI0;			///< 建图时刻惯导系下的原始点
    Vec3T p_L0 = R_LtoI.inverse() * (p_I0 - p_LinI);	///< 建图时刻雷达系下的原始点

    T dist;
    if (measurement_.geo_type == GeometryType::Plane)
	{
	  // 计算点到平面的距离
      Vec3T norm = (measurement_.geo_plane.template cast<T>()).head(3);
      dist = p_L0.dot(norm) + T(measurement_.geo_plane[3]);
    }
	else
	{
      // omit item 1 =: 1.0 / measurement_.geo_normal.norm()
      dist = (p_L0 - measurement_.geo_point.template cast<T>())
                 .cross(measurement_.geo_normal.template cast<T>())
                 .norm();
    }

	// 输出残差
    Eigen::Map<Eigen::Matrix<T, 1, 1>> residuals(sResiduals);
    residuals.template block<1, 1>(0, 0) = Eigen::Matrix<T, 1, 1>(dist);
    residuals = T(weight_) * residuals;

    return true;
  }

 private:
  PointCorrespondence measurement_;		///< 关联点
  SplineMeta<SplineOrder> spline_meta_;	///< 样条节点
  double weight_;						///< 观测值权重
  double inv_dt_;						///< 1.0/dt
};

class PointPlaneFactor {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  PointPlaneFactor(const PointCorrespondence& pc,
                   const SplineMeta<SplineOrder>& spline_meta,
                   double lidar_weight, bool opt_laser_param = false)
      : measurement_(pc),
        spline_meta_(spline_meta),
        lidar_weight_(lidar_weight),
        opt_laser_param_(opt_laser_param)
  {
    inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
  }

  template <class T>
  bool operator()(T const* const* sKnots, T* sResiduals) const {
    using SO3T = Sophus::SO3<T>;
    using Vec2T = Eigen::Matrix<T, 2, 1>;
    using Vec3T = Eigen::Matrix<T, 3, 1>;
    using Vec6T = Eigen::Matrix<T, 6, 1>;

    /// 确定参数位置
    int PARAM_R_LtoI = 0;		// IMU/LiD姿态
    int PARAM_p_LinI = 1;		// IMU/LiD位置
    int PARAM_t_offset = 2;		// LiD外参
    int PARAM_laser_param = 3;  // LiD内参6
    size_t Kont_offset = 2 * spline_meta_.NumParameters();	// 位姿和姿态索引

	// 雷达时延改正
    T t_offset = sKnots[Kont_offset + PARAM_t_offset][0];
    T t = T(measurement_.t_point) + t_offset;

	// 计算插值节点的索引
	T u;			 // 归一化时间
    size_t R_offset; // 旋转开始节点位置
    size_t P_offset; // 位置开始节点位置
    spline_meta_.ComputeSplineIndex(t, R_offset, u);
    P_offset = R_offset + spline_meta_.NumParameters();

	// 内插当前时刻的IMU姿态
    SO3T R_IkToG;
    CeresSplineHelperJet<T, SplineOrder>::template evaluate_lie<Sophus::SO3>(
        sKnots + R_offset, u, inv_dt_, &R_IkToG);

	// 内插当前时刻的IMU位置
    Vec3T p_IkinG;
    CeresSplineHelperJet<T, SplineOrder>::template evaluate<3, 0>(
        sKnots + P_offset, u, inv_dt_, &p_IkinG);

	// 外参
    Eigen::Map<const SO3T> R_LtoI(sKnots[Kont_offset + PARAM_R_LtoI]);
    Eigen::Map<const Vec3T> p_LinI(sKnots[Kont_offset + PARAM_p_LinI]);

	// 计算雷达的位姿
    SO3T R_LkToG = R_IkToG * R_LtoI;
    Vec3T p_LkinG = R_IkToG * p_LinI + p_IkinG;

	///< 雷达系下的原始点
    Vec3T p_Lk;
    if (opt_laser_param_) {
      Vec3T point_raw = measurement_.point_raw.template cast<T>();
      std::vector<T> laser_param;
      laser_param.resize(6);
      for (int i = 0; i < 6; ++i) {
        laser_param.at(i) = sKnots[Kont_offset + PARAM_laser_param + i][0];
      }
      LiDARIntrinsic::GetCalibrated<T>(laser_param, point_raw, p_Lk);
      // if (u > 0.5 && u < 0.50002) {
      //  std::cout << measurement_.point.transpose() << "; "
      //            << p_Lk.transpose() << std::endl;
      //}
    } 
	else
	{
      p_Lk = measurement_.point.template cast<T>();
    }

	// 全局系下的点
    Vec3T p_Gk = R_LkToG * p_Lk + p_LkinG;

	// 计算点到平面的距离
    Vec3T norm = (measurement_.geo_plane.template cast<T>()).head(3);
    T distance = p_Gk.dot(norm) + T(measurement_.geo_plane[3]);

	// 输出残差
    Eigen::Map<Eigen::Matrix<T, 1, 1>> residuals(sResiduals);
    residuals.template block<1, 1>(0, 0) = Eigen::Matrix<T, 1, 1>(distance);
    residuals = T(lidar_weight_) * residuals;

    return true;
  }

 private:
  PointCorrespondence measurement_;		///< 关联点
  SplineMeta<SplineOrder> spline_meta_;	///< 样条节点
  double lidar_weight_;					///< 观测值权重
  double inv_dt_;						///< 1.0/dt
  bool opt_laser_param_;				///< 是否优化雷达参数
};


// 雷达先验位姿因子
class LiDARPoseFactor {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  LiDARPoseFactor(const PoseData& pose_data,
                  const SplineMeta<SplineOrder>& spline_meta, double pos_weight,
                  double rot_weight)
      : pose_data_(pose_data), spline_meta_(spline_meta), pos_weight_(pos_weight), rot_weight_(rot_weight) {
    inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
  }

  template <class T>
  bool operator()(T const* const* sKnots, T* sResiduals) const {
    using Vec3T = Eigen::Matrix<T, 3, 1>;
    using Vec6T = Eigen::Matrix<T, 6, 1>;
    using SO3T = Sophus::SO3<T>;
    using Tangent = typename Sophus::SO3<T>::Tangent;

	// sKnots R_ItoW/p_IinW/R_LtoI/p_LinI/time_offset
    size_t Kont_offset = 2 * spline_meta_.NumParameters();
    T t_offset = sKnots[Kont_offset + 2][0];
    T t = T(pose_data_.timestamp) + t_offset;

    size_t R_offset;  // should be zero if not estimate time offset
    size_t P_offset;
    T u;
    spline_meta_.ComputeSplineIndex(t, R_offset, u);
    P_offset = R_offset + spline_meta_.NumParameters();

	// IMU姿态
    SO3T R_IkToG;
    CeresSplineHelperJet<T, SplineOrder>::template evaluate_lie<Sophus::SO3>(sKnots + R_offset, u, inv_dt_, &R_IkToG);

	// IMU位置
    Vec3T p_IkinG;
    CeresSplineHelperJet<T, SplineOrder>::template evaluate<3, 0>(sKnots + P_offset, u, inv_dt_, &p_IkinG);

	// 雷达外参
    Eigen::Map<SO3T const> const R_LtoI(sKnots[Kont_offset]);
    Eigen::Map<Vec3T const> const p_LinI(sKnots[Kont_offset + 1]);

	// 雷达位姿
    SO3T R_LkToG = R_IkToG * R_LtoI;
    Vec3T p_LkinG = R_IkToG * p_LinI + p_IkinG;

    Eigen::Map<Vec6T> residuals(sResiduals);
    residuals.template block<3, 1>(0, 0) = T(rot_weight_) * (R_LkToG * pose_data_.orientation.inverse()).log();
    residuals.template block<3, 1>(3, 0) = T(pos_weight_) * (p_LkinG - pose_data_.position);
    return true;
  }

 private:
  PoseData pose_data_;					///< 雷达位姿
  SplineMeta<SplineOrder> spline_meta_;	///< 样条节点
  double pos_weight_;					///< 位置权重
  double rot_weight_;					///< 姿态权重
  double inv_dt_;						///< 1.0/dt
};

///< 雷达先验姿态约束
class LiDAROrientationFactor {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  LiDAROrientationFactor(const PoseData& pose_data,
                         const SplineMeta<SplineOrder>& spline_meta,
                         double rot_weight)
      : pose_data_(pose_data), spline_meta_(spline_meta), rot_weight_(rot_weight) {
    inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
  }

  template <class T>
  bool operator()(T const* const* sKnots, T* sResiduals) const {
    using Vec3T = Eigen::Matrix<T, 3, 1>;
    using SO3T = Sophus::SO3<T>;
    using Tangent = typename Sophus::SO3<T>::Tangent;

    T t = T(pose_data_.timestamp);

    size_t R_offset;  // should be zero if not estimate time offset
    T u;
    spline_meta_.ComputeSplineIndex(t, R_offset, u);

	// 计算IMU的姿态 sKnots R_ItoW/R_LtoI
    SO3T R_IkToG;
    CeresSplineHelperJet<T, SplineOrder>::template evaluate_lie<Sophus::SO3>(sKnots + R_offset, u, inv_dt_, &R_IkToG);

	// 雷达外参
    int Kont_offset = spline_meta_.NumParameters();
    Eigen::Map<SO3T const> const R_LtoI(sKnots[Kont_offset]); // 雷达外参

    SO3T R_LkToG = R_IkToG * R_LtoI;
    Eigen::Map<Tangent> residuals(sResiduals);
    residuals = T(rot_weight_) * (R_LkToG * pose_data_.orientation.inverse()).log();
    return true;
  }

 private:
  PoseData pose_data_;					///< 雷达位姿
  SplineMeta<SplineOrder> spline_meta_; ///< 样条节点
  double rot_weight_;					///< 观测值权重
  double inv_dt_;						///< 1.0/dt
};

class SO3KnotFactor {
 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  SO3KnotFactor(const Sophus::SO3d& so3_knots, double so3_knot_weight)
      : so3_knots_(so3_knots), so3_knot_weight_(so3_knot_weight) {}

  template <class T>
  bool operator()(T const* const* sKnots, T* sResiduals) const {
    using Vec3T = Eigen::Matrix<T, 3, 1>;
    using SO3T = Sophus::SO3<T>;
    using Tangent = typename Sophus::SO3<T>::Tangent;

	// sKnots R_ItoW/R_LtoI
    Eigen::Map<SO3T const> const knot(sKnots[0]);   // 雷达姿态
    Eigen::Map<SO3T const> const R_LtoI(sKnots[1]); // 雷达外参

    Eigen::Map<Tangent> residuals(sResiduals);
    residuals = T(so3_knot_weight_) * ((knot * R_LtoI * so3_knots_.inverse()).log());

    return true;
  }

 private:
  Sophus::SO3d so3_knots_;
  double so3_knot_weight_;
};

class PointLocalToGlobalFactor {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	PointLocalToGlobalFactor(const lpostk::PointObs obs, const SplineMeta<SplineOrder>& spline_meta, double weight)
		: obs_(obs), spline_meta_(spline_meta), weight_(weight)
	{
		inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
		num_residuals_ = 3;
	}

	// sKnots 姿态节点 位置节点 外参R 外参T 陆标点
	template <class T>
	bool operator()(T const* const* sKnots, T* sResiduals) const {
		using SO3T = Sophus::SO3<T>;
		using Vec2T = Eigen::Matrix<T, 2, 1>;
		using Vec3T = Eigen::Matrix<T, 3, 1>;
		using Vec6T = Eigen::Matrix<T, 6, 1>;

		if (num_residuals_ != 1 && num_residuals_ != 3)
			return false;

		/// 确定参数位置
		int PARAM_R_LtoI = 0;		// 外参R
		int PARAM_p_LinI = 1;		// 外参T
		int PARAM_MapPoint = 2;		// 陆标点
		size_t Kont_offset = 2 * spline_meta_.NumParameters();	// 位姿和姿态节点个数

		// 当前时间
		T t = T(obs_.timestamp);

		// 计算插值节点的索引
		size_t R_offset;
		size_t P_offset;
		T u;
		spline_meta_.ComputeSplineIndex(t, R_offset, u);	// 姿态
		P_offset = R_offset + spline_meta_.NumParameters(); // 位置

		// 内插当前时刻的IMU姿态
		SO3T R_IkToG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate_lie<Sophus::SO3>(sKnots + R_offset, u, inv_dt_, &R_IkToG);

		// 内插当前时刻的IMU位置
		Vec3T p_IkinG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate<3, 0>(sKnots + P_offset, u, inv_dt_, &p_IkinG);

		Eigen::Map<const SO3T>	R_LtoI(sKnots[Kont_offset + PARAM_R_LtoI]); // 外参R
		Eigen::Map<const Vec3T> p_LinI(sKnots[Kont_offset + PARAM_p_LinI]); // 外参T
		Eigen::Map<const Vec3T> MapPoint(sKnots[Kont_offset + PARAM_MapPoint]); // 陆标点

		// 计算雷达的位姿
		SO3T R_LkToG = R_IkToG * R_LtoI;
		Vec3T p_LkinG = R_IkToG * p_LinI + p_IkinG;

		///< 雷达系下的原始点
		Vec3T p_Lk = Vec3T(T(obs_.XYZ[0]), T(obs_.XYZ[1]), T(obs_.XYZ[2]));

		// 全局系下的点
		Vec3T p_Gk = R_LkToG * p_Lk + p_LkinG;

		// 点与点距离
		T	dist = (p_Gk - MapPoint).norm();

		///////////////////////////////////输出残差//////////////////////////////////////////////
		if (num_residuals_ == 1)
		{
			sResiduals[0] = T(weight_) * (p_Gk - MapPoint).norm();
			return true;
		}
		else if (num_residuals_ == 3)
		{
			Eigen::Map<Vec3T> residuals(sResiduals);
			residuals.template block<3, 1>(0, 0) = T(weight_) * (p_Gk - MapPoint);
			return true;
		}
		///////////////////////////////////输出残差//////////////////////////////////////////////

		return false;
	}

	void set_num_residuals(const int num_residuals)
	{
		num_residuals_ = num_residuals;
	}

	int  num_residuals()
	{
		return num_residuals_;
	}

private:
	lpostk::PointObs	obs_;
	double	weight_;						///< 观测值权重
	double	inv_dt_;						///< 1.0/dt
	SplineMeta<SplineOrder> spline_meta_;	///< 样条节点
	int			num_residuals_;				///< 残差块大小
};

class PointLocalToLocalFactor {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	PointLocalToLocalFactor(const lpostk::PointObs obs1, const lpostk::PointObs obs2, const SplineMeta<SplineOrder>& spline_meta, double weight)
		: obs1_(obs1), obs2_(obs2), spline_meta_(spline_meta), weight_(weight)
	{
		inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
		num_residuals_ = 6;
	}

	// sKnots 姿态节点 位置节点 外参R 外参T 陆标点
	template <class T>
	bool operator()(T const* const* sKnots, T* sResiduals) const {
		using SO3T = Sophus::SO3<T>;
		using Vec2T = Eigen::Matrix<T, 2, 1>;
		using Vec3T = Eigen::Matrix<T, 3, 1>;
		using Vec4T = Eigen::Matrix<T, 4, 1>;
		using Vec6T = Eigen::Matrix<T, 6, 1>;

		if (num_residuals_ != 2 && num_residuals_ != 6)
			return false;

		// 确定参数位置
		int PARAM_R_LtoI = 0;		// 外参R
		int PARAM_p_LinI = 1;		// 外参T
		size_t Kont_offset = 2 * spline_meta_.NumParameters();	// 位姿和姿态节点个数

		/////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////第一帧位姿/////////////////////////////////
		T t1 = T(obs1_.timestamp);
		size_t R1_offset, P1_offset;
		T u1;
		spline_meta_.ComputeSplineIndex(t1, R1_offset, u1);	 // 姿态
		P1_offset = R1_offset + spline_meta_.NumParameters(); // 位置

		SO3T R1_IkToG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate_lie<Sophus::SO3>(sKnots + R1_offset, u1, inv_dt_, &R1_IkToG);
		//std::cout << std::left << R1_IkToG.matrix() << std::endl;

		Vec3T p1_IkinG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate<3, 0>(sKnots + P1_offset, u1, inv_dt_, &p1_IkinG);
		//std::cout << std::left << p1_IkinG << std::endl;
		//////////////////////////////////////第一帧位姿/////////////////////////////////

		//////////////////////////////////////第二帧位姿/////////////////////////////////
		T t2 = T(obs2_.timestamp);
		size_t R2_offset, P2_offset;
		T u2;
		spline_meta_.ComputeSplineIndex(t2, R2_offset, u2);   // 姿态
		P2_offset = R2_offset + spline_meta_.NumParameters(); // 位置

		SO3T R2_IkToG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate_lie<Sophus::SO3>(sKnots + R2_offset, u2, inv_dt_, &R2_IkToG);
		//std::cout << std::left << R2_IkToG.matrix() << std::endl;

		Vec3T p2_IkinG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate<3, 0>(sKnots + P2_offset, u2, inv_dt_, &p2_IkinG);
		//std::cout << std::left << p2_IkinG << std::endl;
		//////////////////////////////////////第二帧位姿/////////////////////////////////

		Eigen::Map<const SO3T>	R_LtoI(sKnots[Kont_offset + PARAM_R_LtoI]); // 外参R
		Eigen::Map<const Vec3T> p_LinI(sKnots[Kont_offset + PARAM_p_LinI]); // 外参T
		//std::cout << std::left << R_LtoI.matrix() << std::endl;
		//std::cout << std::left << p_LinI << std::endl;

		// 相对位姿
		SO3T R_I1ToI2 = R2_IkToG.inverse() * R1_IkToG;
		Vec3T p_I1inI2 = R2_IkToG.inverse() * (p1_IkinG - p2_IkinG);
		SO3T R_I2ToI1 = R1_IkToG.inverse() * R2_IkToG;
		Vec3T p_I2inI1 = R1_IkToG.inverse() * (p2_IkinG - p1_IkinG);

		////////////////////////////////////转换点//////////////////////////////////////////////
		// 第一帧的起点转到第二帧
		Vec3T p1_L1_L1 = Vec3T(T(obs1_.XYZ[0]), T(obs1_.XYZ[1]), T(obs1_.XYZ[2])); ///< 第一帧雷达系下的原始点
		Vec3T p1_L1_I1 = R_LtoI * p1_L1_L1 + p_LinI;				///< 第一帧惯导系下的原始点
		Vec3T p1_L1_I2 = R_I1ToI2 * p1_L1_I1 + p_I1inI2;			///< 第二帧惯导系下的原始点
		Vec3T p1_L1_L2 = R_LtoI.inverse() * (p1_L1_I2 - p_LinI);	///< 第二帧雷达系下的原始点

		// 第二帧的起点转到第一帧
		Vec3T p1_L2_L2 = Vec3T(T(obs2_.XYZ[0]), T(obs2_.XYZ[1]), T(obs2_.XYZ[2])); ///< 第二帧雷达系下的原始点
		Vec3T p1_L2_I2 = R_LtoI * p1_L2_L2 + p_LinI;				///< 第二帧惯导系下的原始点
		Vec3T p1_L2_I1 = R_I2ToI1 * p1_L2_I2 + p_I2inI1;			///< 第一帧惯导系下的原始点
		Vec3T p1_L2_L1 = R_LtoI.inverse() * (p1_L2_I1 - p_LinI);	///< 第一帧雷达系下的原始点
		////////////////////////////////////转换点//////////////////////////////////////////////

		///////////////////////////////////输出残差//////////////////////////////////////////////
		if (num_residuals_ == 2)
		{
			Eigen::Map<Vec2T> residuals(sResiduals);
			residuals(0) = T(weight_) * (p1_L1_L2 - p1_L2_L2).norm();
			residuals(1) = T(weight_) * (p1_L2_L1 - p1_L1_L1).norm();
			return true;
		}
		else if (num_residuals_ == 6)
		{
			Eigen::Map<Vec6T> residuals(sResiduals);
			residuals.template block<3, 1>(0, 0) = T(weight_) * (p1_L1_L2 - p1_L2_L2);
			residuals.template block<3, 1>(3, 0) = T(weight_) * (p1_L2_L1 - p1_L1_L1);
			return true;
		}
		///////////////////////////////////输出残差//////////////////////////////////////////////

		return false;
	}

	void set_num_residuals(const int num_residuals)
	{
		num_residuals_ = num_residuals;
	}

	int  num_residuals()
	{
		return num_residuals_;
	}

private:
	lpostk::PointObs	obs1_;
	lpostk::PointObs	obs2_;
	double	weight_;						///< 观测值权重
	double	inv_dt_;						///< 1.0/dt
	SplineMeta<SplineOrder> spline_meta_;	///< 样条节点
	int			num_residuals_;				///< 残差块大小
};


class LineLocalToGlobalFactor {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	LineLocalToGlobalFactor(const lpostk::LineObs obs, const SplineMeta<SplineOrder>& spline_meta, double dst_weight, const double ang_weight = 1.0 / (5.0 * M_PI / 180.0))
		: obs_(obs), spline_meta_(spline_meta), dst_weight_(dst_weight), ang_weight_(ang_weight)
	{
		inv_dt_ = 1.0 / spline_meta.segments.begin()->dt;
	}

	// sKnots 姿态节点 位置节点 外参R 外参T 全局直线参数
	template <class T>
	bool operator()(T const* const* sKnots, T* sResiduals) const {
		using SO3T = Sophus::SO3<T>;
		using Vec2T = Eigen::Matrix<T, 2, 1>;
		using Vec3T = Eigen::Matrix<T, 3, 1>;
		using Vec6T = Eigen::Matrix<T, 6, 1>;
		using Mat3T = Eigen::Matrix<T, 3, 3>;

		////////////////////////////////////////参数位置///////////////////////////////////
		int PARAM_R_LtoI = 0;		// 外参R
		int PARAM_p_LinI = 1;		// 外参T
		int PARAM_Alpha = 2;		// 直线尺度参数
		int PARAM_Rotation = 3;		// 直线旋转参数
		size_t Kont_offset = 2 * spline_meta_.NumParameters();	// 位姿和姿态节点个数
		////////////////////////////////////////参数位置///////////////////////////////////
		
		//////////////////////////////////////当前帧IMU位姿/////////////////////////////////
		T t = T(obs_.timestamp);
		size_t R_offset, P_offset;
		T u;
		spline_meta_.ComputeSplineIndex(t, R_offset, u);	 // 姿态
		P_offset = R_offset + spline_meta_.NumParameters(); // 位置

		SO3T R_IkToG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate_lie<Sophus::SO3>(sKnots + R_offset, u, inv_dt_, &R_IkToG);

		Vec3T p_IkinG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate<3, 0>(sKnots + P_offset, u, inv_dt_, &p_IkinG);
		//////////////////////////////////////当前帧IMU位姿/////////////////////////////////

		////////////////////////////////////////雷达外参///////////////////////////////////
		Eigen::Map<const SO3T>	R_LtoI(sKnots[Kont_offset + PARAM_R_LtoI]); // 外参R
		Eigen::Map<const Vec3T> p_LinI(sKnots[Kont_offset + PARAM_p_LinI]); // 外参T
		////////////////////////////////////////雷达外参///////////////////////////////////

		////////////////////////////////////////直线参数///////////////////////////////////
		T alpha = sKnots[Kont_offset + PARAM_Alpha][0];	// 尺度参数
		Eigen::Map<const SO3T>	R_Line(sKnots[Kont_offset + PARAM_Rotation]); // 旋转参数
		Mat3T Rotation = R_Line.matrix();				// 旋转矩阵
		Vec3T N0 = Rotation.block<3, 1>(0, 0);			// 单位向量
		Vec3T P0 = alpha * Rotation.block<3, 1>(0, 1);	// 最近点
		////////////////////////////////////////直线参数///////////////////////////////////

		//////////////////////////////////////当前帧LiD位姿/////////////////////////////////
		SO3T R_LkToG = R_IkToG * R_LtoI;
		Vec3T p_LkinG = R_IkToG * p_LinI + p_IkinG;
		//////////////////////////////////////当前帧LiD位姿/////////////////////////////////

		////////////////////////////////////转换点//////////////////////////////////////////////
		Vec3T p1_Lk = Vec3T(T(obs_.P1[0]), T(obs_.P1[1]), T(obs_.P1[2])); // 起点
		Vec3T p1_Gk = R_LkToG * p1_Lk + p_LkinG;
		Vec3T p2_Lk = Vec3T(T(obs_.P2[0]), T(obs_.P2[1]), T(obs_.P2[2])); // 终点
		Vec3T p2_Gk = R_LkToG * p2_Lk + p_LkinG;
		////////////////////////////////////转换点//////////////////////////////////////////////

		///////////////////////////////////计算距离//////////////////////////////////////////////
		T	d1 = (p1_Gk - P0).cross(N0).norm();
		T	d2 = (p2_Gk - P0).cross(N0).norm();
		///////////////////////////////////计算距离//////////////////////////////////////////////

		///////////////////////////////////输出残差//////////////////////////////////////////////
		if (num_residuals_ == 2)
		{
			Eigen::Map<Vec2T> residuals(sResiduals);
			//residuals.template block<1, 1>(0, 0) = T(dst_weight_) * d1;
			//residuals.template block<1, 1>(1, 0) = T(dst_weight_) * d2;
			return true;
		}
		///////////////////////////////////输出残差//////////////////////////////////////////////

		return false;
	}

	void set_num_residuals(const int num_residuals)
	{
		num_residuals_ = num_residuals;
	}

	int  num_residuals()
	{
		return num_residuals_;
	}

private:
	lpostk::LineObs		obs_;						///< 平面特征
	double		dst_weight_;				///< 距离权重
	double		ang_weight_;				///< 角度权重
	double		inv_dt_;					///< 1.0/dt
	SplineMeta<SplineOrder> spline_meta_;	///< 样条节点
	int			num_residuals_;				///< 残差块大小
};

class LineLocalToLocalFactor {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	LineLocalToLocalFactor(const lpostk::LineObs obs1, const lpostk::LineObs obs2, 
		const SplineMeta<SplineOrder>& spline_meta, double dst_weight, const double ang_weight)
		: obs1_(obs1), obs2_(obs2), spline_meta_(spline_meta), dst_weight_(dst_weight), ang_weight_(ang_weight)
	{
		inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
	}

	// sKnots 姿态节点 位置节点 外参R 外参T
	template <class T>
	bool operator()(T const* const* sKnots, T* sResiduals) const {
		using SO3T = Sophus::SO3<T>;
		using Vec2T = Eigen::Matrix<T, 2, 1>;
		using Vec3T = Eigen::Matrix<T, 3, 1>;
		using Vec4T = Eigen::Matrix<T, 4, 1>;
		using Vec6T = Eigen::Matrix<T, 6, 1>;

		// 确定参数位置
		int PARAM_R_LtoI = 0;		// 外参R
		int PARAM_p_LinI = 1;		// 外参T
		size_t Kont_offset = 2 * spline_meta_.NumParameters();	// 位姿和姿态节点个数

		/////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////第一帧位姿/////////////////////////////////
		T t1 = T(obs1_.timestamp);
		size_t R1_offset, P1_offset;
		T u1;
		spline_meta_.ComputeSplineIndex(t1, R1_offset, u1);	 // 姿态
		P1_offset = R1_offset + spline_meta_.NumParameters(); // 位置

		SO3T R1_IkToG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate_lie<Sophus::SO3>(sKnots + R1_offset, u1, inv_dt_, &R1_IkToG);

		Vec3T p1_IkinG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate<3, 0>(sKnots + P1_offset, u1, inv_dt_, &p1_IkinG);
		//////////////////////////////////////第一帧位姿/////////////////////////////////

		//////////////////////////////////////第二帧位姿/////////////////////////////////
		T t2 = T(obs2_.timestamp);
		size_t R2_offset, P2_offset;
		T u2;
		spline_meta_.ComputeSplineIndex(t2, R2_offset, u2);   // 姿态
		P2_offset = R2_offset + spline_meta_.NumParameters(); // 位置

		SO3T R2_IkToG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate_lie<Sophus::SO3>(sKnots + R2_offset, u2, inv_dt_, &R2_IkToG);

		Vec3T p2_IkinG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate<3, 0>(sKnots + P2_offset, u2, inv_dt_, &p2_IkinG);
		//////////////////////////////////////第二帧位姿/////////////////////////////////

		Eigen::Map<const SO3T>	R_LtoI(sKnots[Kont_offset + PARAM_R_LtoI]); // 外参R
		Eigen::Map<const Vec3T> p_LinI(sKnots[Kont_offset + PARAM_p_LinI]); // 外参T

		// 相对位姿
		SO3T R_I1ToI2 = R2_IkToG.inverse() * R1_IkToG;
		Vec3T p_I1inI2 = R2_IkToG.inverse() * (p1_IkinG - p2_IkinG);
		SO3T R_I2ToI1 = R1_IkToG.inverse() * R2_IkToG;
		Vec3T p_I2inI1 = R1_IkToG.inverse() * (p2_IkinG - p1_IkinG);

		////////////////////////////////////转换点//////////////////////////////////////////////
		// 第一帧的起点转到第二帧
		Vec3T p1_L1_L1 = Vec3T(T(obs1_.P1[0]), T(obs1_.P1[1]), T(obs1_.P1[2])); ///< 第一帧雷达系下的原始点
		Vec3T p1_L1_I1 = R_LtoI * p1_L1_L1 + p_LinI;				///< 第一帧惯导系下的原始点
		Vec3T p1_L1_I2 = R_I1ToI2 * p1_L1_I1 + p_I1inI2;			///< 第二帧惯导系下的原始点
		Vec3T p1_L1_L2 = R_LtoI.inverse() * (p1_L1_I2 - p_LinI);	///< 第二帧雷达系下的原始点

		// 第一帧的终点转到第二帧
		Vec3T p2_L1_L1 = Vec3T(T(obs1_.P1[0]), T(obs1_.P1[1]), T(obs1_.P1[2])); ///< 第一帧雷达系下的原始点
		Vec3T p2_L1_I1 = R_LtoI * p2_L1_L1 + p_LinI;				///< 第一帧惯导系下的原始点
		Vec3T p2_L1_I2 = R_I1ToI2 * p2_L1_I1 + p_I1inI2;			///< 第二帧惯导系下的原始点
		Vec3T p2_L1_L2 = R_LtoI.inverse() * (p2_L1_I2 - p_LinI);	///< 第二帧雷达系下的原始点

		// 第二帧的起点转到第一帧
		Vec3T p1_L2_L2 = Vec3T(T(obs2_.P1[0]), T(obs2_.P1[1]), T(obs2_.P1[2])); ///< 第二帧雷达系下的原始点
		Vec3T p1_L2_I2 = R_LtoI * p1_L2_L2 + p_LinI;				///< 第二帧惯导系下的原始点
		Vec3T p1_L2_I1 = R_I2ToI1 * p1_L2_I2 + p_I2inI1;			///< 第一帧惯导系下的原始点
		Vec3T p1_L2_L1 = R_LtoI.inverse() * (p1_L2_I1 - p_LinI);	///< 第一帧雷达系下的原始点

		// 第二帧的终点转到第一帧
		Vec3T p2_L2_L2 = Vec3T(T(obs2_.P2[0]), T(obs2_.P2[1]), T(obs2_.P2[2])); ///< 第二帧雷达系下的原始点
		Vec3T p2_L2_I2 = R_LtoI * p2_L2_L2 + p_LinI;				///< 第二帧惯导系下的原始点
		Vec3T p2_L2_I1 = R_I2ToI1 * p2_L2_I2 + p_I2inI1;			///< 第一帧惯导系下的原始点
		Vec3T p2_L2_L1 = R_LtoI.inverse() * (p2_L2_I1 - p_LinI);	///< 第一帧雷达系下的原始点
		////////////////////////////////////转换点//////////////////////////////////////////////

		///////////////////////////////////计算距离//////////////////////////////////////////////
		Vec3T p1_p2_L1 = p2_L1_L1 - p1_L1_L1;
		Vec3T p1_p2_L2 = p2_L2_L2 - p1_L2_L2;
		T	d1_L1 = (p1_L2_L1 - p1_L1_L1).cross(p1_p2_L1).norm() / p1_p2_L1.norm();
		T	d2_L1 = (p2_L2_L1 - p1_L1_L1).cross(p1_p2_L1).norm() / p1_p2_L1.norm();
		T	d1_L2 = (p1_L1_L2 - p1_L2_L2).cross(p1_p2_L2).norm() / p1_p2_L2.norm();
		T	d2_L2 = (p2_L1_L2 - p1_L2_L2).cross(p1_p2_L2).norm() / p1_p2_L2.norm();
		///////////////////////////////////计算距离//////////////////////////////////////////////

		///////////////////////////////////输出残差//////////////////////////////////////////////
		if (num_residuals_ == 2)
		{
			sResiduals[0] = T(dst_weight_) * d1_L1;
			sResiduals[1] = T(dst_weight_) * d2_L1;
			return true;
		}
		else if (num_residuals_ == 4)
		{
			sResiduals[0] = T(dst_weight_) * d1_L1;
			sResiduals[1] = T(dst_weight_) * d2_L1;
			sResiduals[2] = T(dst_weight_) * d1_L2;
			sResiduals[3] = T(dst_weight_) * d2_L2;
			return true;
		}
		///////////////////////////////////输出残差//////////////////////////////////////////////

		return false;
	}


	void set_num_residuals(const int num_residuals)
	{
		num_residuals_ = num_residuals;
	}

	int  num_residuals()
	{
		return num_residuals_;
	}

private:
	lpostk::LineObs		obs1_;				///< 直线特征
	lpostk::LineObs		obs2_;				///< 直线特征
	double		dst_weight_;				///< 距离权重
	double		ang_weight_;				///< 角度权重
	double		inv_dt_;					///< 1.0/dt
	SplineMeta<SplineOrder> spline_meta_;	///< 样条节点
	int			num_residuals_;				///< 残差块大小
};

class PoleLocalToGlobalFactor {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	PoleLocalToGlobalFactor(const lpostk::PoleObs obs, const SplineMeta<SplineOrder>& spline_meta, double dst_weight, const double ang_weight = 1.0)
		: obs_(obs), spline_meta_(spline_meta), dst_weight_(dst_weight), ang_weight_(ang_weight)
	{
		inv_dt_ = 1.0 / spline_meta.segments.begin()->dt;
		num_residuals_ = 2;
	}

	// sKnots 姿态节点 位置节点 外参R 外参T
	template <class T>
	bool operator()(T const* const* sKnots, T* sResiduals) const {
		using SO3T = Sophus::SO3<T>;
		using Vec2T = Eigen::Matrix<T, 2, 1>;
		using Vec3T = Eigen::Matrix<T, 3, 1>;
		using Vec6T = Eigen::Matrix<T, 6, 1>;
		using Mat3T = Eigen::Matrix<T, 3, 3>;

		////////////////////////////////////////参数位置///////////////////////////////////
		int PARAM_R_LtoI = 0;		// 外参R
		int PARAM_p_LinI = 1;		// 外参T
		int PARAM_Alpha = 2;		// 轴线尺度参数
		int PARAM_Rotation = 3;		// 轴线旋转参数
		size_t Kont_offset = 2 * spline_meta_.NumParameters();	// 位姿和姿态节点个数
		////////////////////////////////////////参数位置///////////////////////////////////

		//////////////////////////////////////当前帧IMU位姿/////////////////////////////////
		T t = T(obs_.timestamp);
		size_t R_offset, P_offset;
		T u;
		spline_meta_.ComputeSplineIndex(t, R_offset, u);	 // 姿态
		P_offset = R_offset + spline_meta_.NumParameters(); // 位置

		SO3T R_IkToG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate_lie<Sophus::SO3>(sKnots + R_offset, u, inv_dt_, &R_IkToG);

		Vec3T p_IkinG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate<3, 0>(sKnots + P_offset, u, inv_dt_, &p_IkinG);
		//////////////////////////////////////当前帧IMU位姿/////////////////////////////////

		////////////////////////////////////////雷达外参///////////////////////////////////
		Eigen::Map<const SO3T>	R_LtoI(sKnots[Kont_offset + PARAM_R_LtoI]); // 外参R
		Eigen::Map<const Vec3T> p_LinI(sKnots[Kont_offset + PARAM_p_LinI]); // 外参T
		////////////////////////////////////////雷达外参///////////////////////////////////

		////////////////////////////////////////直线参数///////////////////////////////////
		T alpha = sKnots[Kont_offset + PARAM_Alpha][0]; // 尺度参数
		Eigen::Map<const SO3T>	R_Line(sKnots[Kont_offset + PARAM_Rotation]); // 旋转参数
		Mat3T Rotation = R_Line.matrix();
		Vec3T N0 = Rotation.block<3, 1>(0, 0);			// 单位向量
		Vec3T P0 = alpha * Rotation.block<3, 1>(0, 1);	// 最近点
		////////////////////////////////////////直线参数///////////////////////////////////

		//////////////////////////////////////当前帧LiD位姿/////////////////////////////////
		SO3T R_LkToG = R_IkToG * R_LtoI;
		Vec3T p_LkinG = R_IkToG * p_LinI + p_IkinG;
		//////////////////////////////////////当前帧LiD位姿/////////////////////////////////

		////////////////////////////////////转换点//////////////////////////////////////////////
		// 当前帧的起点转到全局系
		Vec3T p1_Lk = Vec3T(T(obs_.P1[0]), T(obs_.P1[1]), T(obs_.P1[2]));
		Vec3T p1_Gk = R_LkToG * p1_Lk + p_LkinG;				///< 全局系

		// 当前帧的终点转到第一帧
		Vec3T p2_Lk = Vec3T(T(obs_.P2[0]), T(obs_.P2[1]), T(obs_.P2[2]));
		Vec3T p2_Gk = R_LkToG * p2_Lk + p_LkinG;				///< 第二帧惯导系下的原始点
		////////////////////////////////////转换点//////////////////////////////////////////////

		///////////////////////////////////计算距离//////////////////////////////////////////////
		T	d1 = (p1_Gk - P0).cross(N0).norm();
		T	d2 = (p2_Gk - P0).cross(N0).norm();
		///////////////////////////////////计算距离//////////////////////////////////////////////
		
		///////////////////////////////////输出残差//////////////////////////////////////////////
		if (num_residuals_ == 2)
		{
			Eigen::Map<Vec2T> residuals(sResiduals);
			//residuals.template block<1, 1>(0, 0) = T(dst_weight_) * d1;
			//residuals.template block<1, 1>(1, 0) = T(dst_weight_) * d2;
			return true;
		}
		///////////////////////////////////输出残差//////////////////////////////////////////////

		return false;
	}

	void set_num_residuals(const int num_residuals)
	{
		num_residuals_ = num_residuals;
	}

	int  num_residuals()
	{
		return num_residuals_;
	}

private:
	lpostk::PoleObs		obs_;						///< 平面特征
	double		dst_weight_;				///< 距离权重
	double		ang_weight_;				///< 角度权重
	double		inv_dt_;					///< 1.0/dt
	SplineMeta<SplineOrder> spline_meta_;	///< 样条节点
	int			num_residuals_;				///< 残差块大小
};

class PoleLocalToLocalFactor {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	PoleLocalToLocalFactor(const lpostk::PoleObs obs1, const lpostk::PoleObs obs2, 
			const SplineMeta<SplineOrder>& spline_meta, double dst_weight, const double ang_weight)
		: obs1_(obs1), obs2_(obs2), spline_meta_(spline_meta), dst_weight_(dst_weight), ang_weight_(ang_weight)
	{
		inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
		num_residuals_ = 2;
	}

	// sKnots 姿态节点 位置节点 外参R 外参T 陆标点
	template <class T>
	bool operator()(T const* const* sKnots, T* sResiduals) const {
		using SO3T = Sophus::SO3<T>;
		using Vec2T = Eigen::Matrix<T, 2, 1>;
		using Vec3T = Eigen::Matrix<T, 3, 1>;
		using Vec4T = Eigen::Matrix<T, 4, 1>;
		using Vec6T = Eigen::Matrix<T, 6, 1>;

		if (num_residuals_ != 2 && num_residuals_ != 4)
			return false;

		// 确定参数位置
		int PARAM_R_LtoI = 0;		// 外参R
		int PARAM_p_LinI = 1;		// 外参T
		size_t Kont_offset = 2 * spline_meta_.NumParameters();	// 位姿和姿态节点个数

		/////////////////////////////////////////////////////////////////////////////////
		//////////////////////////////////////第一帧位姿/////////////////////////////////
		T t1 = T(obs1_.timestamp);
		size_t R1_offset, P1_offset;
		T u1;
		spline_meta_.ComputeSplineIndex(t1, R1_offset, u1);	 // 姿态
		P1_offset = R1_offset + spline_meta_.NumParameters(); // 位置

		SO3T R1_IkToG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate_lie<Sophus::SO3>(sKnots + R1_offset, u1, inv_dt_, &R1_IkToG);

		Vec3T p1_IkinG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate<3, 0>(sKnots + P1_offset, u1, inv_dt_, &p1_IkinG);
		//////////////////////////////////////第一帧位姿/////////////////////////////////

		//////////////////////////////////////第二帧位姿/////////////////////////////////
		T t2 = T(obs2_.timestamp);
		size_t R2_offset, P2_offset;
		T u2;
		spline_meta_.ComputeSplineIndex(t2, R2_offset, u2);   // 姿态
		P2_offset = R2_offset + spline_meta_.NumParameters(); // 位置

		SO3T R2_IkToG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate_lie<Sophus::SO3>(sKnots + R1_offset, u2, inv_dt_, &R2_IkToG);

		Vec3T p2_IkinG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate<3, 0>(sKnots + P1_offset, u2, inv_dt_, &p2_IkinG);
		//////////////////////////////////////第二帧位姿/////////////////////////////////

		Eigen::Map<const SO3T>	R_LtoI(sKnots[Kont_offset + PARAM_R_LtoI]); // 外参R
		Eigen::Map<const Vec3T> p_LinI(sKnots[Kont_offset + PARAM_p_LinI]); // 外参T

		// 相对位姿
		SO3T R_I1ToI2 = R2_IkToG.inverse() * R1_IkToG;
		Vec3T p_I1inI2 = R2_IkToG.inverse() * (p1_IkinG - p2_IkinG);
		SO3T R_I2ToI1 = R1_IkToG.inverse() * R2_IkToG;
		Vec3T p_I2inI1 = R1_IkToG.inverse() * (p2_IkinG - p1_IkinG);

		////////////////////////////////////转换点//////////////////////////////////////////////
		// 第一帧的起点转到第二帧
		Vec3T p1_L1_L1 = Vec3T(T(obs1_.P1[0]), T(obs1_.P1[1]), T(obs1_.P1[2])); ///< 第一帧雷达系下的原始点
		Vec3T p1_L1_I1 = R_LtoI * p1_L1_L1 + p_LinI;				///< 第一帧惯导系下的原始点
		Vec3T p1_L1_I2 = R_I1ToI2 * p1_L1_I1 + p_I1inI2;			///< 第二帧惯导系下的原始点
		Vec3T p1_L1_L2 = R_LtoI.inverse() * (p1_L1_I2 - p_LinI);	///< 第二帧雷达系下的原始点

		// 第一帧的终点转到第二帧
		Vec3T p2_L1_L1 = Vec3T(T(obs1_.P1[0]), T(obs1_.P1[1]), T(obs1_.P1[2])); ///< 第一帧雷达系下的原始点
		Vec3T p2_L1_I1 = R_LtoI * p2_L1_L1 + p_LinI;				///< 第一帧惯导系下的原始点
		Vec3T p2_L1_I2 = R_I1ToI2 * p2_L1_I1 + p_I1inI2;			///< 第二帧惯导系下的原始点
		Vec3T p2_L1_L2 = R_LtoI.inverse() * (p2_L1_I2 - p_LinI);	///< 第二帧雷达系下的原始点

		// 第二帧的起点转到第一帧
		Vec3T p1_L2_L2 = Vec3T(T(obs2_.P1[0]), T(obs2_.P1[1]), T(obs2_.P1[2])); ///< 第二帧雷达系下的原始点
		Vec3T p1_L2_I2 = R_LtoI * p1_L2_L2 + p_LinI;				///< 第二帧惯导系下的原始点
		Vec3T p1_L2_I1 = R_I2ToI1 * p1_L2_I2 + p_I2inI1;			///< 第一帧惯导系下的原始点
		Vec3T p1_L2_L1 = R_LtoI.inverse() * (p1_L2_I1 - p_LinI);	///< 第一帧雷达系下的原始点

		// 第二帧的终点转到第一帧
		Vec3T p2_L2_L2 = Vec3T(T(obs2_.P2[0]), T(obs2_.P2[1]), T(obs2_.P2[2])); ///< 第二帧雷达系下的原始点
		Vec3T p2_L2_I2 = R_LtoI * p2_L2_L2 + p_LinI;				///< 第二帧惯导系下的原始点
		Vec3T p2_L2_I1 = R_I2ToI1 * p2_L2_I2 + p_I2inI1;			///< 第一帧惯导系下的原始点
		Vec3T p2_L2_L1 = R_LtoI.inverse() * (p2_L2_I1 - p_LinI);	///< 第一帧雷达系下的原始点
		////////////////////////////////////转换点//////////////////////////////////////////////

		///////////////////////////////////计算距离//////////////////////////////////////////////
		Vec3T p1_p2_L1 = p2_L1_L1 - p1_L1_L1;
		Vec3T p1_p2_L2 = p2_L2_L2 - p1_L2_L2;
		T	d1_L1 = ((p1_L2_L1 - p1_L1_L1).cross(p1_p2_L1)).norm() / p1_p2_L1.norm();
		T	d2_L1 = ((p2_L2_L1 - p1_L1_L1).cross(p1_p2_L1)).norm() / p1_p2_L1.norm();
		T	d1_L2 = ((p1_L1_L2 - p1_L2_L2).cross(p1_p2_L2)).norm() / p1_p2_L2.norm();
		T	d2_L2 = ((p2_L1_L2 - p1_L2_L2).cross(p1_p2_L2)).norm() / p1_p2_L2.norm();
		///////////////////////////////////计算距离//////////////////////////////////////////////

		///////////////////////////////////输出残差//////////////////////////////////////////////
		if (num_residuals_ == 2)
		{
			// 单向距离
			sResiduals[0] = T(dst_weight_) * d1_L1;
			sResiduals[1] = T(dst_weight_) * d2_L1;
			return true;
		}
		else if (num_residuals_ == 4)
		{
			// 双向距离
			sResiduals[0] = T(dst_weight_) * d1_L1;
			sResiduals[1] = T(dst_weight_) * d2_L1;
			sResiduals[2] = T(dst_weight_) * d1_L2;
			sResiduals[3] = T(dst_weight_) * d2_L2;
			return true;
		}
		///////////////////////////////////输出残差//////////////////////////////////////////////

		return false;
	}

	void set_num_residuals(const int num_residuals)
	{
		num_residuals_ = num_residuals;
	}

	int  num_residuals()
	{
		return num_residuals_;
	}

private:
	lpostk::PoleObs		obs1_;				///< 柱状特征
	lpostk::PoleObs		obs2_;				///< 柱状特征
	double		dst_weight_;				///< 距离权重
	double		ang_weight_;				///< 角度权重
	double		inv_dt_;					///< 1.0/dt
	SplineMeta<SplineOrder> spline_meta_;	///< 样条节点
	int			num_residuals_;				///< 残差块大小
};

// 当前帧平面全局参数
class PlaneLocalToGlobalCPFactor {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		PlaneLocalToGlobalCPFactor(const lpostk::PlaneObs obs, const SplineMeta<SplineOrder>& spline_meta, double dst_weight, const double ang_weight)
		: obs_(obs), spline_meta_(spline_meta), dst_weight_(dst_weight), ang_weight_(ang_weight)
	{
		inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
	}

	// sKnots 姿态节点 位置节点 外参R 外参T 陆标点
	template <class T>
	bool operator()(T const* const* sKnots, T* sResiduals) const {
		using SO3T = Sophus::SO3<T>;
		using Vec2T = Eigen::Matrix<T, 2, 1>;
		using Vec3T = Eigen::Matrix<T, 3, 1>;
		using Vec6T = Eigen::Matrix<T, 6, 1>;

		/// 确定参数位置
		int PARAM_R_LtoI = 0;		// 外参R
		int PARAM_p_LinI = 1;		// 外参T
		int PARAM_Plane = 2;		// 陆标点
		size_t Kont_offset = 2 * spline_meta_.NumParameters();	// 位姿和姿态节点个数

		// 当前时间
		T t = T(obs_.timestamp);

		// 计算插值节点的索引
		size_t R_offset;
		size_t P_offset;
		T u;
		spline_meta_.ComputeSplineIndex(t, R_offset, u);	// 姿态
		P_offset = R_offset + spline_meta_.NumParameters(); // 位置

		// 内插当前时刻的IMU姿态
		SO3T R_IkToG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate_lie<Sophus::SO3>(
			sKnots + R_offset, u, inv_dt_, &R_IkToG);

		// 内插当前时刻的IMU位置
		Vec3T p_IkinG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate<3, 0>(
			sKnots + P_offset, u, inv_dt_, &p_IkinG);

		Eigen::Map<const SO3T>	R_LtoI(sKnots[Kont_offset + PARAM_R_LtoI]); // 外参R
		Eigen::Map<const Vec3T> p_LinI(sKnots[Kont_offset + PARAM_p_LinI]); // 外参T
		Eigen::Map<const Vec3T> CPW(sKnots[Kont_offset + PARAM_Plane]);	// 全局系下的平面最近点
		Vec3T NW = CPW / CPW.norm(); // 全局系下的单位法向量

		// 计算雷达的位姿
		SO3T R_LkToG = R_IkToG * R_LtoI;
		Vec3T p_LkinG = R_IkToG * p_LinI + p_IkinG;

		// 平面中心点
		Vec3T p_Lk = Vec3T(T(obs_.Center[0]), T(obs_.Center[1]), T(obs_.Center[2]));
		Vec3T p_Gk = R_LkToG * p_Lk + p_LkinG;

		// 平面法向量
		Vec3T N_Lk = Vec3T(T(obs_.Equation[0]), T(obs_.Equation[1]), T(obs_.Equation[2]));
		Vec3T N_Gk = R_LkToG * N_Lk;

		// 点到平面的距离
		T	distance = p_Gk.dot(CPW) / CPW.norm() - CPW.norm();

		// 法向量夹角
		T	theta = acos(NW.dot(N_Gk)); // 这个角度在0~90读附近或者为90~180之间,如果为后者要调整为0~90
		if (theta > M_PI_2)	theta = M_PI - acos(NW.dot(N_Gk));

		///////////////////////////////////输出残差//////////////////////////////////////////////
		if (num_residuals_ == 1)
		{
			sResiduals[0] = T(dst_weight_) * distance;
		}
		else if (num_residuals_ == 2)
		{
			sResiduals[0] = T(dst_weight_) * distance;
			sResiduals[1] = T(ang_weight_) * theta;
		}
		///////////////////////////////////输出残差//////////////////////////////////////////////

		return true;
	}

	void set_num_residuals(const int num_residuals)
	{
		num_residuals_ = num_residuals;
	}

	int  num_residuals()
	{
		return num_residuals_;
	}

private:
	lpostk::PlaneObs	obs_;						///< 平面特征
	double		dst_weight_;				///< 距离权重
	double		ang_weight_;				///< 角度权重
	double		inv_dt_;					///< 1.0/dt
	SplineMeta<SplineOrder> spline_meta_;	///< 样条节点
	int			num_residuals_;				///< 残差块大小
};

//两帧之间的平面特征
class PlaneLocalToLocalFactor {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	PlaneLocalToLocalFactor(const lpostk::PlaneObs obs1, const lpostk::PlaneObs obs2, 
		const SplineMeta<SplineOrder>& spline_meta, double dst_weight, const double ang_weight)
		: obs1_(obs1), obs2_(obs2), spline_meta_(spline_meta), dst_weight_(dst_weight), ang_weight_(ang_weight)
	{
		inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
		num_residuals_ = 2;
	}

	// sKnots 姿态节点 位置节点 外参R 外参T 陆标点
	template <class T>
	bool operator()(T const* const* sKnots, T* sResiduals) const {
		using SO3T = Sophus::SO3<T>;
		using Vec2T = Eigen::Matrix<T, 2, 1>;
		using Vec3T = Eigen::Matrix<T, 3, 1>;
		using Vec4T = Eigen::Matrix<T, 4, 1>;
		using Vec6T = Eigen::Matrix<T, 6, 1>;

		if (num_residuals_ != 2 && num_residuals_ != 3 && num_residuals_ != 4)
			return false;

		///// 确定参数位置
		int PARAM_R_LtoI = 0;		// 外参R
		int PARAM_p_LinI = 1;		// 外参T
		size_t Kont_offset = 2 * spline_meta_.NumParameters();	// 位姿和姿态节点个数

		////////////////////////////////////第一帧位姿/////////////////////////////////
		T t1 = T(obs1_.timestamp);
		size_t R1_offset, P1_offset;
		T u1;
		spline_meta_.ComputeSplineIndex(t1, R1_offset, u1);	 // 姿态
		P1_offset = R1_offset + spline_meta_.NumParameters(); // 位置

		SO3T R1_IkToG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate_lie<Sophus::SO3>(sKnots + R1_offset, u1, inv_dt_, &R1_IkToG);

		Vec3T p1_IkinG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate<3, 0>(sKnots + P1_offset, u1, inv_dt_, &p1_IkinG);
		////////////////////////////////////第一帧位姿/////////////////////////////////

		////////////////////////////////////第二帧位姿/////////////////////////////////
		T t2 = T(obs2_.timestamp);
		size_t R2_offset, P2_offset;
		T u2;
		spline_meta_.ComputeSplineIndex(t2, R2_offset, u2);   // 姿态
		P2_offset = R2_offset + spline_meta_.NumParameters(); // 位置

		SO3T R2_IkToG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate_lie<Sophus::SO3>(sKnots + R1_offset, u2, inv_dt_, &R2_IkToG);

		Vec3T p2_IkinG;
		CeresSplineHelperJet<T, SplineOrder>::template evaluate<3, 0>(sKnots + P1_offset, u2, inv_dt_, &p2_IkinG);
		////////////////////////////////////第二帧位姿/////////////////////////////////

		Eigen::Map<const SO3T>	R_LtoI(sKnots[Kont_offset + PARAM_R_LtoI]); // 外参R
		Eigen::Map<const Vec3T> p_LinI(sKnots[Kont_offset + PARAM_p_LinI]); // 外参T

		// 第一帧的中心点转到第二帧
		SO3T R_I1ToI2 = R2_IkToG.inverse() * R1_IkToG;
		Vec3T p_I1inI2 = R2_IkToG.inverse() * (p1_IkinG - p2_IkinG);
		Vec3T p_L1_L1 = Vec3T(T(obs1_.Center[0]), T(obs1_.Center[1]), T(obs1_.Center[2])); ///< 第一帧雷达系下的原始点
		Vec3T p_L1_I1 = R_LtoI * p_L1_L1 + p_LinI;				///< 第一帧惯导系下的原始点
		Vec3T p_L1_I2 = R_I1ToI2 * p_L1_I1 + p_I1inI2;			///< 第二帧惯导系下的原始点
		Vec3T p_L1_L2 = R_LtoI.inverse() * (p_L1_I2 - p_LinI);	///< 第二帧雷达系下的原始点

		// 第二帧的中心点转到第一帧
		SO3T R_I2ToI1 = R1_IkToG.inverse() * R2_IkToG;
		Vec3T p_I2inI1 = R1_IkToG.inverse() * (p2_IkinG - p1_IkinG);
		Vec3T p_L2_L2 = Vec3T(T(obs2_.Center[0]), T(obs2_.Center[1]), T(obs2_.Center[2])); ///< 第二帧雷达系下的原始点
		Vec3T p_L2_I2 = R_LtoI * p_L2_L2 + p_LinI;				///< 第二帧惯导系下的原始点
		Vec3T p_L2_I1 = R_I2ToI1 * p_L2_I2 + p_I2inI1;			///< 第一帧惯导系下的原始点
		Vec3T p_L2_L1 = R_LtoI.inverse() * (p_L2_I1 - p_LinI);	///< 第一帧雷达系下的原始点

		// 法向量
		Vec3T norm_L1_L1 = Vec3T(T(obs1_.Equation[0]), T(obs1_.Equation[1]), T(obs1_.Equation[2])); // 单向法向量
		Vec3T norm_L2_L2 = Vec3T(T(obs2_.Equation[0]), T(obs2_.Equation[1]), T(obs2_.Equation[2])); // 单位法向量
		Vec3T norm_L1_L2 = R_I1ToI2 * norm_L1_L1;
		Vec3T norm_L2_L1 = R_I2ToI1 * norm_L2_L2;

		// 点到平面的距离(双向)
		T distance_L1 = p_L2_L1.dot(norm_L1_L1) + T(obs1_.Equation[3]);
		T distance_L2 = p_L1_L2.dot(norm_L2_L2) + T(obs2_.Equation[3]);

		// 法向量夹角
		T	theta_L1 = acos(norm_L1_L1.dot(norm_L2_L1));
		if (theta_L1 > M_PI_2)	theta_L1 = M_PI - acos(norm_L1_L1.dot(norm_L2_L1)); // 调整到0~90
		T	theta_L2 = acos(norm_L2_L2.dot(norm_L1_L2));
		if (theta_L2 > M_PI_2)	theta_L2 = M_PI - acos(norm_L2_L2.dot(norm_L1_L2)); // 调整到0~90

		///////////////////////////////////输出残差//////////////////////////////////////////////
		if (num_residuals_ == 2)
		{
			// 单向距离 + 单向法向量夹角
			Eigen::Map<Vec2T> residuals(sResiduals);
			//residuals.template block<1, 1>(0, 0) = T(dst_weight_) * distance_L1;
			//residuals.template block<1, 1>(1, 0) = T(ang_weight_) * theta_L1;
			return true;
		}
		else if (num_residuals_ == 3)
		{
			// 双向距离 + 单向法向量夹角
			Eigen::Map<Vec3T> residuals(sResiduals);
			//residuals.template block<1, 1>(0, 0) = T(dst_weight_) * distance_L1;
			//residuals.template block<1, 1>(1, 0) = T(dst_weight_) * distance_L2;
			//residuals.template block<1, 1>(2, 0) = T(ang_weight_) * theta_L1;
		}
		else if (num_residuals_ == 4)
		{
			// 双向距离 + 双向法向量夹角
			Eigen::Map<Vec4T> residuals(sResiduals);
			//residuals.template block<1, 1>(0, 0) = T(dst_weight_) * distance_L1;
			//residuals.template block<1, 1>(1, 0) = T(dst_weight_) * distance_L2;
			//residuals.template block<1, 1>(2, 0) = T(ang_weight_) * theta_L1;
			//residuals.template block<1, 1>(3, 0) = T(ang_weight_) * theta_L2;
			return true;
		}
		///////////////////////////////////输出残差//////////////////////////////////////////////

		return true;
	}

	void set_num_residuals(const int num_residuals)
	{
		num_residuals_ = num_residuals;
	}

	int  num_residuals()
	{
		return num_residuals_;
	}

private:
	lpostk::PlaneObs	obs1_;				///< 平面特征
	lpostk::PlaneObs	obs2_;				///< 平面特征
	double		dst_weight_;				///< 距离权重
	double		ang_weight_;				///< 角度权重
	double		inv_dt_;					///< 1.0/dt
	SplineMeta<SplineOrder> spline_meta_;	///< 样条节点
	int			num_residuals_;				///< 残差块大小
};


}  // namespace liso

#endif
