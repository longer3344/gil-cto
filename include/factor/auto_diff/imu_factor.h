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

#ifndef IMU_FACTOR_H
#define IMU_FACTOR_H

#include <basalt/spline/ceres_spline_helper.h>
#include <basalt/spline/spline_segment.h>
#include <ceres/ceres.h>
#include <sensor_data/imu_data.h>
#include <sophus/so3.hpp>

namespace liso {
	using namespace basalt;


	// 加计/陀螺先验零偏约束
	class IMUBiasPriorFactor {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		IMUBiasPriorFactor(const double weight)
		{
			bias_.setZero();
			weight3_(0) = weight3_(1) = weight3_(2) = weight;
		}

		IMUBiasPriorFactor(const Eigen::Vector3d weight3) : weight3_(weight3)
		{
			bias_.setZero();
		}

		IMUBiasPriorFactor(const double bias, const double weight)
		{
			bias_(0) = bias_(1) = bias_(2) = bias;
			weight3_(0) = weight3_(1) = weight3_(2) = weight;
		}

		IMUBiasPriorFactor(const Eigen::Vector3d& bias, const double weight) : bias_(bias)
		{
			weight3_(0) = weight3_(1) = weight3_(2) = weight;
		}

		IMUBiasPriorFactor(const Eigen::Vector3d& bias, const Eigen::Vector3d& weight3) : bias_(bias), weight3_(weight3)
		{
		}

		template <class T>
		bool operator()(T const* const* sKnots, T* sResiduals) const {

			using Vec3T = Eigen::Matrix<T, 3, 1>;

			// step 1 : 获取零偏
			Eigen::Map<const Vec3T> bias(sKnots[0]);

			// step 2 : 计算残差
			Eigen::Map<Vec3T> residuals(sResiduals);
			residuals = bias - bias_.template cast<T>();
			residuals(0) = T(weight3_(0)) * residuals(0);
			residuals(1) = T(weight3_(1)) * residuals(1);
			residuals(2) = T(weight3_(2)) * residuals(2);
			//residuals = T(weight_) * residuals;
			return true;
		}

	private:
		Eigen::Vector3d	bias_;			///< 零偏
		Eigen::Vector3d	weight3_;		///< 权重, 1.0 / std
	};


	// 零偏之间先验约束
	class IMUBiasBetweenFactor {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		IMUBiasBetweenFactor(const double weight)
		{
			weight3_(0) = weight3_(1) = weight3_(2);
		}

		IMUBiasBetweenFactor(const Eigen::Vector3d weight3) : weight3_(weight3)
		{

		}

		template <class T>
		bool operator()(T const* const* sKnots, T* sResiduals) const {

			using Vec3T = Eigen::Matrix<T, 3, 1>;

			// step 1 : 获取零偏
			Eigen::Map<const Vec3T> bias1(sKnots[0]); // 零偏1
			Eigen::Map<const Vec3T> bias2(sKnots[1]); // 零偏2

			// step 2 : 计算残差
			Eigen::Map<Vec3T> residuals(sResiduals);
			residuals = bias2 - bias1;
			residuals(0) = T(weight3_(0)) * residuals(0);
			residuals(1) = T(weight3_(1)) * residuals(1);
			residuals(2) = T(weight3_(2)) * residuals(2);
			return true;
		}

	private:
		Eigen::Vector3d weight3_;		///< 权重
	};

	// 陀螺数据因子
	class GyroscopeWithConstantBiasFactor : public CeresSplineHelper<SplineOrder> {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		GyroscopeWithConstantBiasFactor(const IMUData& imu_data, const SplineMeta<SplineOrder>& spline_meta, double weight)
			: imu_data_(imu_data), spline_meta_(spline_meta), weight_(weight) {
			inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
		}

		template <class T>
		bool operator()(T const* const* sKnots, T* sResiduals) const {
			using Vec3T = Eigen::Matrix<T, 3, 1>;
			using Tangent = typename Sophus::SO3<T>::Tangent;

			// step 1 : 计算开始节点 sKnots R_ItoW/bias
			size_t R_offset;
			double u;
			spline_meta_.ComputeSplineIndex(imu_data_.timestamp, R_offset, u);

			// step 2 : 计算角速度
			Tangent rot_vel;
			CeresSplineHelper<SplineOrder>::template evaluate_lie<T, Sophus::SO3>(sKnots + R_offset, u, inv_dt_, nullptr, &rot_vel, nullptr);

			// step 3 : 计算残差
			size_t Kont_offset = spline_meta_.NumParameters();
			Eigen::Map<const Vec3T> gyro_bias(sKnots[Kont_offset]);		// 陀螺零偏

			Eigen::Map<Vec3T> residuals(sResiduals);
			residuals = rot_vel - imu_data_.gyro.template cast<T>() + gyro_bias;
			residuals = T(weight_) * residuals;
			return true;
		}

	private:
		IMUData imu_data_; ///< IMU数据
		SplineMeta<SplineOrder> spline_meta_; ///< 样条节点
		double weight_;    ///< 陀螺数据权重
		double inv_dt_;    ///< 样条时间间隔
	};


	// 加计数据因子
	class AccelerometerWithConstantBiasFactor : public CeresSplineHelper<SplineOrder> {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		AccelerometerWithConstantBiasFactor(const IMUData& imu_data, const SplineMeta<SplineOrder>& spline_meta, double weight,
			const bool blocal, const Eigen::Vector3d gravity = Eigen::Vector3d(0.0, 0.0, 0.0))
			: imu_data_(imu_data), spline_meta_(spline_meta), weight_(weight), blocal_(blocal), gravity_(gravity) {
			inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
		}

		template <class T>
		bool operator()(T const* const* sKnots, T* sResiduals) const {
			using Vec3T = Eigen::Matrix<T, 3, 1>;
			using SO3T = Sophus::SO3<T>;
			using Tangent = typename Sophus::SO3<T>::Tangent;

			// step 1 : 计算当前时刻节点位置
			// sKnots R_ItoW/p_IinW/acce_bias
			size_t R_offset;
			size_t P_offset;
			double u;
			spline_meta_.ComputeSplineIndex(imu_data_.timestamp, R_offset, u);
			P_offset = R_offset + spline_meta_.NumParameters();

			// step 2 : 计算当前时刻位姿
			// step 2.1 : 计算姿态
			SO3T R_w_i; // Rbe
			CeresSplineHelper<SplineOrder>::template evaluate_lie<T, Sophus::SO3>(sKnots + R_offset, u, inv_dt_, &R_w_i, nullptr); // 姿态

			// step 2.2 : 计算位置和加速度
			Vec3T position;
			Vec3T velocity;
			Vec3T accel_w;
			CeresSplineHelper<SplineOrder>::template evaluate<T, 3, 0>(sKnots + P_offset, u, inv_dt_, &position);
			CeresSplineHelper<SplineOrder>::template evaluate<T, 3, 1>(sKnots + P_offset, u, inv_dt_, &velocity);
			CeresSplineHelper<SplineOrder>::template evaluate<T, 3, 2>(sKnots + P_offset, u, inv_dt_, &accel_w);

			size_t Kont_offset = 2 * spline_meta_.NumParameters();
			Eigen::Map<const Vec3T> acce_bias(sKnots[Kont_offset]);	// 加计零偏

			// 导出加速度
			Vec3T wie = Eigen::Vector3d(0.0, 0.0, 7.2921151467E-5).template cast<T>(); // 地球自转
			Vec3T coriolis = T(2.0)*wie.cross(velocity); // 科室力
			if (gravity_.norm() == 0.0)
			{
				fprintf(stdout, "Error in AccelerometerWithConstantBiasFactor gravity not set !!");
				return false;
			}

			//Vec3T a_b = R_w_i.inverse() * (accel_w - gravity_);
			Vec3T a_b = R_w_i.inverse() * (accel_w + coriolis - gravity_);


			// 加计残差
			Vec3T acce_residuals = a_b - imu_data_.accel.template cast<T>() + acce_bias;
			//std::cout << "acce_residuals " << std::endl << acce_residuals << std::endl;

			Eigen::Map<Vec3T> residuals(sResiduals);
			residuals.template block<3, 1>(0, 0) = T(weight_) * acce_residuals;
			return true;
		}

	private:
		IMUData imu_data_; ///< IMU数据
		SplineMeta<SplineOrder> spline_meta_; ///< 样条节点
		double weight_;    ///< 陀螺数据权重
		double inv_dt_;    ///< 样条时间间隔
		bool   blocal_;
		Eigen::Vector3d gravity_; // 当前时刻重力
	};

	// 陀螺和加计因子
	class GyroAcceWithConstantBiasFactor : public CeresSplineHelper<SplineOrder> {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		GyroAcceWithConstantBiasFactor(const IMUData& imu_data, const SplineMeta<SplineOrder>& spline_meta,
				const double gyro_weight, const double acce_weight,
				const bool blocal, const Eigen::Vector3d gravity = Eigen::Vector3d(0.0, 0.0, 0.0))
			: imu_data_(imu_data), spline_meta_(spline_meta),
			gyro_weight_(gyro_weight), acce_weight_(acce_weight),
			blocal_(blocal), gravity_(gravity) {
			inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
		}

		template <class T>
		bool operator()(T const* const* sKnots, T* sResiduals) const
		{
			using Vec2T = Eigen::Matrix<T, 2, 1>;
			using Vec3T = Eigen::Matrix<T, 3, 1>;
			using Vec6T = Eigen::Matrix<T, 6, 1>;
			using Vec9T = Eigen::Matrix<T, 9, 1>;
			using Mat3T = Eigen::Matrix<T, 3, 3>;
			using SO3T = Sophus::SO3<T>;
			using Tangent = typename Sophus::SO3<T>::Tangent;

			// 地球自转角速度
			Eigen::Vector3d Wie(0.0, 0.0, 7.2921151467E-5);

			// 局部坐标系下的轨迹要检查重力加速度          
			if (!blocal_ && gravity_.norm() <= 1.0)
				return false;

			static double timestamp = 0.0;
			double  aaa = std::remainder(imu_data_.timestamp, 1.0);
			if (0 && timestamp != imu_data_.timestamp && fabs(fmod(imu_data_.timestamp + 1.0e-06, 1.0)) < 0.001)
				printf("imu_data_ %10.3lf, %lf\n", imu_data_.timestamp, aaa);
			timestamp = imu_data_.timestamp;
			//printf("imu_data_ %14.6lf\n", imu_data_.timestamp);

			// step 1 : 计算当前时刻节点位置
			// sKnots R_ItoW/p_IinW/gyro_bias/acce_bias/g_refine/Mw/Ma/Aw/q_WtoA
			size_t R_offset;  // should be zero if not estimate time offset
			size_t P_offset;
			double u;
			spline_meta_.ComputeSplineIndex(imu_data_.timestamp, R_offset, u);
			P_offset = R_offset + spline_meta_.NumParameters();

			// step 2 : 计算当前时刻位姿
			// step 2.1 : 计算姿态和角速度
			SO3T R_w_i; // Rbe
			Tangent rot_vel; // 
			CeresSplineHelper<SplineOrder>::template evaluate_lie<T, Sophus::SO3>(sKnots + R_offset, u, inv_dt_, &R_w_i, &rot_vel);

			// step 2.2 : 计算位置和加速度
			Vec3T position;
			Vec3T velocity;
			Vec3T accel_w;
			CeresSplineHelper<SplineOrder>::template evaluate<T, 3, 0>(sKnots + P_offset, u, inv_dt_, &position);
			CeresSplineHelper<SplineOrder>::template evaluate<T, 3, 1>(sKnots + P_offset, u, inv_dt_, &velocity);
			CeresSplineHelper<SplineOrder>::template evaluate<T, 3, 2>(sKnots + P_offset, u, inv_dt_, &accel_w);

			size_t Kont_offset = 2 * spline_meta_.NumParameters();
			Eigen::Map<const Vec3T> gyro_bias(sKnots[Kont_offset]);		// 陀螺零偏
			Eigen::Map<const Vec3T> acce_bias(sKnots[Kont_offset + 1]);	// 加计零偏
			Eigen::Map<const Vec2T> g_refine(sKnots[Kont_offset + 2]);	// 重力方向, 局部坐标系下的俯仰角和横滚角

			auto Mw_vec = sKnots[Kont_offset + 3];     // 陀螺尺度因子(0-2)和非正交误差(3-5)
			auto Ma_vec = sKnots[Kont_offset + 4];     // 加计尺度因子(0-2)和非正交误差(3-5)
			auto Aw_vec = sKnots[Kont_offset + 5];     // 加计对陀螺的影响
			auto q_WtoA_vec = sKnots[Kont_offset + 6]; // 陀螺相对于加计的安置角
			Eigen::Map<SO3T const> const S_WtoA(q_WtoA_vec);
			//fprintf(stdout, "q_WtoA %10.6lf %10.6lf %10.6lf %10.6lf\n", q_WtoA_vec[0], q_WtoA_vec[1], q_WtoA_vec[2], q_WtoA_vec[3]);

			// 陀螺尺度因子和非正交误差
			Mat3T Mw = Mat3T::Zero();
			Mw.diagonal() = Eigen::Map<const Vec3T>(Mw_vec, 3);	// 陀螺尺度因子
			Mw(0, 1) = *(Mw_vec + 3);
			Mw(0, 2) = *(Mw_vec + 4);
			Mw(1, 2) = *(Mw_vec + 5);
			//std::cout << Mw << std::endl << std::endl;

			// 加计尺度因子和非正交误差
			Mat3T Ma = Mat3T::Zero();
			Ma.diagonal() = Eigen::Map<const Vec3T>(Ma_vec, 3);	// 加计尺度因子
			Ma(0, 1) = *(Ma_vec + 3);
			Ma(0, 2) = *(Ma_vec + 4);
			Ma(1, 2) = *(Ma_vec + 5);
			//std::cout << Ma << std::endl << std::endl;

			// 加计对陀螺的影响
			Mat3T Aw = Mat3T::Zero();
			Aw.col(0) = Eigen::Map<const Vec3T>(Aw_vec, 3);
			Aw.col(1) = Eigen::Map<const Vec3T>(Aw_vec + 3, 3);
			Aw.col(2) = Eigen::Map<const Vec3T>(Aw_vec + 6, 3);
			//std::cout << Aw << std::endl << std::endl;

			// 重力
			Vec3T gravity = gravity_factor::refined_gravity(g_refine);
			Vec3T a_b = R_w_i.inverse() * (accel_w + gravity);
			
			Vec3T wie = Eigen::Vector3d(0.0, 0.0, 7.2921151467E-5).template cast<T>();
			Vec3T coriolis = T(2.0)*wie.cross(velocity);
			if (1 && !blocal_) // 全局坐标系下, 扣掉重力加速度
			{
				gravity = gravity_.template cast<T>();
				a_b = R_w_i.inverse() * (accel_w + coriolis - gravity);
			}

			// 计算IMU的角速度和加速度
			Vec3T w_b = (Mw * (S_WtoA * rot_vel) + Aw * a_b).eval();
			a_b = (Ma * a_b).eval(); // 尺度因子和非正交误差

			// 陀螺残差
			Vec3T gyro_residuals = w_b - imu_data_.gyro.template cast<T>() + gyro_bias;
			if (1 && !blocal_)
			{
				//gyro_residuals = w_b + R_w_i.inverse() * Wie - imu_data_.gyro.template cast<T>() + gyro_bias;
			}

			// 加计残差
			Vec3T acce_residuals = a_b - imu_data_.accel.template cast<T>() + acce_bias;  

			//std::cout << "gyro_residuals " << std::endl << gyro_residuals << std::endl;
			//std::cout << "acce_residuals " << std::endl << acce_residuals << std::endl;

			Eigen::Map<Vec6T> residuals(sResiduals);
			residuals.template block<3, 1>(0, 0) = T(gyro_weight_) * gyro_residuals;
			residuals.template block<3, 1>(3, 0) = T(acce_weight_) * acce_residuals;
			return true;
		}

	private:
		IMUData imu_data_;
		SplineMeta<SplineOrder> spline_meta_;
		double gyro_weight_;
		double acce_weight_;
		double inv_dt_;
		bool   blocal_;
		Eigen::Vector3d gravity_;
	};

	// IMU位姿因子
	class IMUPoseFactor : public CeresSplineHelper<SplineOrder> {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		IMUPoseFactor(const PoseData& pose_data, const SplineMeta<SplineOrder>& spline_meta,
				double pos_weight, double rot_weight)
			: pose_data_(pose_data), spline_meta_(spline_meta),
			pos_weight_(pos_weight), rot_weight_(rot_weight) {
			inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
		}

		template <class T>
		bool operator()(T const* const* sKnots, T* sResiduals) const {
			using Vec3T = Eigen::Matrix<T, 3, 1>;
			using Vec6T = Eigen::Matrix<T, 6, 1>;
			using SO3T = Sophus::SO3<T>;
			using Tangent = typename Sophus::SO3<T>::Tangent;

			// step 1 : 计算当前时间的样条节点位置
			size_t R_offset; // 位置
			size_t P_offset; // 姿态
			double u;
			spline_meta_.ComputeSplineIndex(pose_data_.timestamp, R_offset, u);
			P_offset = R_offset + spline_meta_.NumParameters();

			// step 2.1: 计算当前时间的姿态
			SO3T R_IkToG;
			CeresSplineHelper<SplineOrder>::template evaluate_lie<T, Sophus::SO3>(sKnots + R_offset, u, inv_dt_, &R_IkToG);

			// step 2.1: 计算当前时间的位置
			Vec3T p_IkinG;
			CeresSplineHelper<SplineOrder>::template evaluate<T, 3, 0>(sKnots + P_offset, u, inv_dt_, &p_IkinG);

			// step 3 : 计算残差
			Eigen::Map<Vec6T> residuals(sResiduals);
			residuals.template block<3, 1>(0, 0) = T(rot_weight_) * (R_IkToG * pose_data_.orientation.inverse()).log();
			residuals.template block<3, 1>(3, 0) = T(pos_weight_) * (p_IkinG - pose_data_.position);
			return true;
		}

	private:
		PoseData pose_data_;
		SplineMeta<SplineOrder> spline_meta_;
		double pos_weight_;
		double rot_weight_;
		double inv_dt_;
	};

	// IMU位置因子
	class IMUPositionFactor : public CeresSplineHelper<SplineOrder> {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		// pose_data - IMU位姿 spline_meta - 当前位姿的样条支撑节点 pos_weight - 位置权重
		IMUPositionFactor(const PoseData& pose_data, const SplineMeta<SplineOrder>& spline_meta, double pos_weight)
			: timestamp_(pose_data.timestamp), position_(pose_data.position), spline_meta_(spline_meta), pos_weight_(pos_weight) {
			inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
		}

		// timestamp - 当前时间 position - IMU位置 spline_meta - 当前位姿的样条支撑节点 pos_weight - 位置权重
		IMUPositionFactor(const double timestamp, const Eigen::Vector3d& position,
			const SplineMeta<SplineOrder>& spline_meta, double pos_weight)
			: timestamp_(timestamp), position_(position), spline_meta_(spline_meta), pos_weight_(pos_weight) {
			inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
		}

		template <class T>
		bool operator()(T const* const* sKnots, T* sResiduals) const {
			using Vec3T = Eigen::Matrix<T, 3, 1>;

			// step 1 : 计算当前时间的样条节点位置
			size_t P_offset;
			double u;
			spline_meta_.ComputeSplineIndex(timestamp_, P_offset, u);

			// step 2 : 计算当前时间的位置
			Vec3T p_IkinG;
			CeresSplineHelper<SplineOrder>::template evaluate<T, 3, 0>(sKnots + P_offset, u, inv_dt_, &p_IkinG);

			// step 3 : 计算残差
			Eigen::Map<Vec3T> residuals(sResiduals);
			residuals = T(pos_weight_) * (p_IkinG - position_);

			return true;
		}

	private:
		double timestamp_;
		Eigen::Vector3d position_;
		SplineMeta<SplineOrder> spline_meta_;
		double pos_weight_;
		double inv_dt_;
	};

	// IMU速度因子
	class IMUVelocityFactor : public CeresSplineHelper<SplineOrder> {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		// pose_data - IMU位姿 spline_meta - 当前位姿的样条支撑节点 vel_weight - 速度权重
		IMUVelocityFactor(const PoseData& pose_data, const SplineMeta<SplineOrder>& spline_meta, double vel_weight)
			: timestamp_(pose_data.timestamp), velocity_(pose_data.velocity), spline_meta_(spline_meta), vel_weight_(vel_weight) {
			inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
		}

		// timestamp - 当前时间 velocity - IMU速度 spline_meta - 当前位姿的样条支撑节点 vel_weight - 速度权重
		IMUVelocityFactor(const double timestamp, const Eigen::Vector3d& velocity,
			const SplineMeta<SplineOrder>& spline_meta, double vel_weight)
			: timestamp_(timestamp), velocity_(velocity), spline_meta_(spline_meta), vel_weight_(vel_weight) {
			inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
		}

		template <class T>
		bool operator()(T const* const* sKnots, T* sResiduals) const {
			using Vec3T = Eigen::Matrix<T, 3, 1>;

			// step 1 : 计算当前时间的样条节点位置
			size_t P_offset;
			double u;
			spline_meta_.ComputeSplineIndex(timestamp_, P_offset, u);

			// step 2 : 计算当前时间的速度
			Vec3T velocity;
			CeresSplineHelper<SplineOrder>::template evaluate<T, 3, 1>(sKnots + P_offset, u, inv_dt_, &velocity);

			// step 3 : 计算残差
			Eigen::Map<Vec3T> residuals(sResiduals);
			residuals = T(vel_weight_) * (velocity - velocity_);

			return true;
		}

	private:
		double timestamp_;
		Eigen::Vector3d velocity_;
		SplineMeta<SplineOrder> spline_meta_;
		double vel_weight_;
		double inv_dt_;
	};

	// IMU姿态因子
	class IMUOrientationFactor : public CeresSplineHelper<SplineOrder> {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

		// pose_data - IMU位姿 spline_meta - 当前位姿的样条支撑节点 rot_weight - 姿态权重
		IMUOrientationFactor(const PoseData& pose_data, const SplineMeta<SplineOrder>& spline_meta, double rot_weight)
			: timestamp_(pose_data.timestamp), orientation_(pose_data.orientation), spline_meta_(spline_meta), rot_weight_(rot_weight) {
			inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
		}

		// timestamp - 当前时间 orientation - IMU姿态 spline_meta - 当前位姿的样条支撑节点 rot_weight - 姿态权重
		IMUOrientationFactor(const double timestamp, const SO3d& orientation,
			const SplineMeta<SplineOrder>& spline_meta, double rot_weight)
			: timestamp_(timestamp), orientation_(orientation), spline_meta_(spline_meta), rot_weight_(rot_weight) {
			inv_dt_ = 1.0 / spline_meta_.segments.begin()->dt;
		}

		template <class T>
		bool operator()(T const* const* sKnots, T* sResiduals) const {
			using SO3T = Sophus::SO3<T>;
			using Tangent = typename Sophus::SO3<T>::Tangent;

			// step 1 : 计算当前时间的样条节点位置
			size_t R_offset;
			double u;
			spline_meta_.ComputeSplineIndex(timestamp_, R_offset, u);

			// step 2 : 计算当前时间的速度
			SO3T R_IkToG;
			CeresSplineHelper<SplineOrder>::template evaluate_lie<T, Sophus::SO3>(sKnots, u, 1, &R_IkToG);

			// step 3 : 计算残差
			Eigen::Map<Tangent> residuals(sResiduals);
			residuals = T(rot_weight_) * ((R_IkToG * orientation_.inverse()).log());

			return true;
		}

	private:
		double timestamp_;
		SO3d orientation_;
		SplineMeta<SplineOrder> spline_meta_;
		double rot_weight_;
		double inv_dt_;
	};

}  // namespace liso

#endif
