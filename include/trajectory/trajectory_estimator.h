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

#ifndef TRAJECTORY_ESTIMATOR_H
#define TRAJECTORY_ESTIMATOR_H

#include <ceres/ceres.h>
#include <ceres/covariance.h>
#include <factor/auto_diff/gnss_factor.h>
#include <factor/auto_diff/imu_factor.h>
#include <factor/auto_diff/lidar_feature_factor.h>
#include <factor/auto_diff/motion_factor.h>
#include <trajectory/se3_trajectory.h>
#include <utils/ceres_callbacks.h>
#include <basalt/spline/ceres_local_param.hpp>
#include <memory>
#include <thread>

 //#include "calib/surfel_association.h"
#include "sensor_data/calibration.h"
#include "sensor_data/imu_data.h"
#include "trajectory/trajectory_estimator_options.h"

#include "LiDARFFStream.h"
using namespace lpostk;

enum LiDMeasType
{
	LOCALTOGLOBAL,
	LOCALTOLOCAL
};

namespace liso {

	class TrajectoryEstimator {
		static ceres::Problem::Options DefaultProblemOptions() {
			ceres::Problem::Options options;
			options.loss_function_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
			options.local_parameterization_ownership = ceres::DO_NOT_TAKE_OWNERSHIP;
			return options;
		}

	public:
		/// @brief constructor
		///
		/// @param trajectory  轨迹节点, 包括位置和姿态
		/// @param calib_param 标定参数
		/// @param option      优化配置参数
		TrajectoryEstimator(std::shared_ptr<Trajectory> trajectory, CalibParamManager::Ptr calib_param, const TrajectoryEstimatorOptions& option)
			: trajectory_(trajectory), calib_param_(calib_param), options_(option) {
			problem_ = std::make_shared<ceres::Problem>(DefaultProblemOptions());
			local_parameterization = new LieLocalParameterization<SO3d>();
			pFrameContainer = NULL;
		}

		void SetProblem(std::shared_ptr<ceres::Problem> problem_in) { problem_ = problem_in; }

		void SetTrajectory(std::shared_ptr<Trajectory> trajectory_in) { trajectory_ = trajectory_in; }

		bool CheckSplineMeta(const SplineMeta<SplineOrder>& spline_meta);

		void AddControlPoints(const SplineMeta<SplineOrder>& spline_meta, std::vector<double*>& vec, bool addPosKont = false);

		void AddLoamMeasurement(const PointCorrespondence& pc, double weight, double huber_loss = 5);

		void AddLiDARPoseMeasurement(const PoseData& pose_data, double rot_weight, double pos_weight);

		void AddLiDARSurfelMeasurement(const PointCorrespondence& pc, double lidar_weight);

		/// @brief 添加陀螺数据
		///
		/// @param[in] imu_data	    IMU数据
		/// @param[in] gyro_weight	陀螺权重
		/// @return void
		void AddIMUGyroMeasurement(const IMUData& imu_data, double gyro_weight);

		/// @brief 添加陀螺数据
		///
		/// @param[in] timestamp	当前时间
		/// @param[in] gyro_data	陀螺数据
		/// @param[in] gyro_weight	陀螺权重
		/// @return void
		void AddIMUGyroMeasurement(const double timestamp, const Eigen::Vector3d& gyro_data, double gyro_weight);

		/// @brief 添加加计数据
		///
		/// @param[in] imu_data	    IMU数据
		/// @param[in] accel_weight	加计权重
		/// @return void
		void AddIMUAccelMeasurement(const IMUData& imu_data, double accel_weight);

		/// @brief 添加加计数据
		///
		/// @param[in] timestamp	当前时间
		/// @param[in] accel_data	加计数据
		/// @param[in] accel_weight	加计权重
		/// @return void
		void AddIMUAccelMeasurement(const double timestamp, const Eigen::Vector3d& accel_data, double accel_weight);


		/// @brief 添加IMU数据
		///
		/// @param[in] imu_data	    IMU数据
		/// @param[in] gyro_weight	陀螺权重
		/// @param[in] accel_weight	加计权重
		/// @return void
		void AddIMUMeasurement(const IMUData& imu_data, double gyro_weight, double accel_weight);

		/// @brief 添加IMU位姿因子
		///
		/// @param[in] pose_data	IMU位姿
		/// @param[in] rot_weight	姿态权重
		/// @param[in] pos_weight	位置权重
		/// @return void
		void AddPoseMeasurement(const PoseData& pose_data, double rot_weight, double pos_weight);

		/// @brief IMU位置
		///
		/// @param[in] pose_data	IMU位姿
		/// @param[in] pos_weight	位置权重
		/// @return void
		void AddPositionMeasurement(const PoseData& pose_data, double pos_weight);

		/// @brief IMU位置
		///
		/// @param[in] timestamp	当前时间
		/// @param[in] position	IMU位置
		/// @param[in] pos_weight	位置权重
		/// @return void
		void AddPositionMeasurement(const double timestamp, const Eigen::Vector3d position, double pos_weight);

		/// @brief IMU速度
		///
		/// @param[in] pose_data	IMU位姿
		/// @param[in] pos_weight	速度权重
		/// @return void
		void AddVelocityMeasurement(const PoseData& pose_data, double vel_weight);

		/// @brief IMU速度
		///
		/// @param[in] timestamp	当前时间
		/// @param[in] Velocity	IMU速度
		/// @param[in] vel_weight	速度权重
		/// @return void
		void AddVelocityMeasurement(const double timestamp, const Eigen::Vector3d Velocity, double vel_weight);

		/// @brief IMU姿态
		///
		/// @param[in] pose_data	IMU位姿
		/// @param[in] pos_weight	姿态权重
		/// @return void
		void AddOrientationMeasurement(const PoseData& pose_data, double rot_weight);

		/// @brief IMU姿态
		///
		/// @param[in] timestamp	 当前时间
		/// @param[in] Orientation IMU姿态
		/// @param[in] rot_weight	 姿态权重
		/// @return void
		void AddOrientationMeasurement(const double timestamp, const SO3d& Orientation, double rot_weight);

		/// @brief 零偏先验约束
		///
		/// @param[in] bias			先验零偏
		/// @param[in] weight		先验零偏的权重
		/// @return void
		void AddAccelBiasPriorConstraint(const Eigen::Vector3d weight);
		void AddAccelBiasPriorConstraint(const double weight);
		void AddAccelBiasPriorConstraint(const Eigen::Vector3d &bias, const Eigen::Vector3d weight);
		void AddAccelBiasPriorConstraint(const Eigen::Vector3d &bias, const double weight);
		void AddAccelBiasPriorConstraint(const double bias, const double weight);
		void AddGyroBiasPriorConstraint(const Eigen::Vector3d weight);
		void AddGyroBiasPriorConstraint(const double weight);
		void AddGyroBiasPriorConstraint(const Eigen::Vector3d &bias, const Eigen::Vector3d weight);
		void AddGyroBiasPriorConstraint(const Eigen::Vector3d &bias, const double weight);
		void AddGyroBiasPriorConstraint(const double bias, const double weight);

		//
		void AddAccelPriorBias(const double biasStd);
		void AddAccelPriorBias(const double bias, const double biasStd);
		void AddGyroPriorBias(const double bias, const double biasStd);
		void AddGyroPriorBias(const double biasStd);

	
		/// @brief 零偏随机游走约束
		///
		/// @param[in] PNSD			零偏随机游走谱密度
		/// @return void
		void AddAccelBiasBetweenConstraint(const Eigen::Vector3d PNSD);
		void AddAccelBiasBetweenConstraint(const double PNSD);
		void AddGyroBiasBetweenConstraint(const Eigen::Vector3d PNSD);
		void AddGyroBiasBetweenConstraint(const double PNSD);
		void AddIMUBiasBetweenConstraint(const Eigen::Vector3d GyrPNSD, const Eigen::Vector3d AccPNSD);
		void AddIMUBiasBetweenConstraint(const double GyrPNSD, const double AccPNSD);

		void AddQuadraticIntegralFactor(double min_time, double max_time, Eigen::Matrix3d weight);

		void AddAngularVelocityConvexHullFactor(double time, double weight);

		void printPointResiduals();
		bool AddPointMeasurements(const LiDMeasType type, const double weight = 20);
		ceres::ResidualBlockId AddPointLocalToGlobalMeasurement(const PointObs &obs, double PW[3], double weight);
		ceres::ResidualBlockId AddPointLocalToLocalMeasurement(const PointObs &obs1, const PointObs &obs2, double weight);

		void printLineResiduals();
		bool AddLineMeasurements(const LiDMeasType type, const double dst_weight, const double ang_weight);
		ceres::ResidualBlockId AddLineLocalToGlobalMeasurement(const LineObs &obs, double *Alpha, liso::SO3d *Rotation, double dst_weight, const double ang_weight);
		ceres::ResidualBlockId AddLineLocalToLocalMeasurement(const LineObs &obs1, const LineObs &obs2, double dst_weight, const double ang_weight);

		void printPoleResiduals();
		bool AddPoleMeasurements(const LiDMeasType type, const double dst_weight, const double ang_weight);
		ceres::ResidualBlockId AddPoleLocalToGlobalMeasurement(const PoleObs &obs, double *Alpha, liso::SO3d *Rotation, double dst_weight, const double ang_weight);
		ceres::ResidualBlockId AddPoleLocalToLocalMeasurement(const PoleObs &obs1, const PoleObs &obs2, double dst_weight, const double ang_weight);

		void printPlaneResiduals();
		bool AddPlaneMeasurements(const LiDMeasType type, const double dst_weight, const double ang_weight);
		ceres::ResidualBlockId AddPlaneLocalToGlobalMeasurement(const PlaneObs &obs, double CPW[3], double dst_weight, const double ang_weight);
		ceres::ResidualBlockId AddPlaneLocalToLocalMeasurement(const PlaneObs &obs1, const PlaneObs &obs2, double dst_weight, const double ang_weight);

		/// @brief GNSS位置/速度约束
		///
		/// @param[in] obs			GNSS解算结果
		/// @param[in] blarm		GNSS杆臂
		/// @param[in] weight		先验位置/速度的权重
		/// @return void
		ceres::ResidualBlockId addGNSSPosMeasurement(const GNSSPVA &obs, const Eigen::Vector3d &blarm, double weight);
		ceres::ResidualBlockId addGNSSVelMeasurement(const GNSSPVA &obs, const Eigen::Vector3d &blarm, double weight);

		void SetTrajectorControlPointVariable(double min_time, double max_time);

		/// Add callback for debug
		void AddCallback(bool needs_state = true, std::string filename = "");

		void ShowParameters();

		ceres::Solver::Summary Solve(int max_iterations = 50, bool progress = true, int num_threads = 1);

		bool getCovariance();

	public:

		/// @brief 锁定轨迹节点
		///
		/// @param[in] lock_P 是否锁定位置节点
		/// @param[in] lock_R 是否锁定姿态节点
		/// @return void
		void LockTrajKnots(bool lock_P, bool lock_R);

		/// @brief 锁定惯导参数
		///
		/// @param[in] lock_ab 是否锁定加计零偏
		/// @param[in] lock_wb 是否锁定陀螺零偏
		/// @param[in] lock_g  是否锁定重力
		/// @return void
		void LockIMUState(bool lock_ab, bool lock_wb, bool lock_g);

		/// @brief 锁定惯导参数
		///
		/// @param[in] lock_P 是否锁定杆臂
		/// @param[in] lock_R 是否锁定安置角
		/// @return void
		void LockExtrinsicParam(bool lock_P, bool lock_R);

		/// @brief 输出分段惯导参数
		///
		/// @return void
		void ShowSegmentParam() const
		{
			for (int i = 0; i < calib_param_->segment_param.size(); i++)
			{
				auto& v = calib_param_->segment_param[i];
				fprintf(stdout, "%4d %13.6lf %13.6lf %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e\n", i, v.starttime, v.endtime,
					v.acce_bias(0), v.acce_bias(1), v.acce_bias(2), v.gyro_bias(0), v.gyro_bias(1), v.gyro_bias(2));
			}
		}

		CalibParamManager::Ptr calib_param_; ///< 标定参数
		TrajectoryEstimatorOptions options_; ///< 优化配置参数
		std::shared_ptr<Trajectory> trajectory_;  ///< 轨迹节点, 包括位置和姿态
		std::shared_ptr<ceres::Problem> problem_; ///< ceres优化问题
		ceres::LocalParameterization* local_parameterization; ///< 旋转矩阵本地参数化
		bool callback_needs_state_; ///< 是否回调
		std::vector<std::unique_ptr<ceres::IterationCallback>> callbacks_; ///< 回调函数
		lpostk::FrameContainer* pFrameContainer; ///< 激光雷达数据
	};

}  // namespace liso

#endif
