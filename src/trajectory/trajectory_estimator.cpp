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

#include <trajectory/trajectory_estimator.h>

namespace liso {

	bool TrajectoryEstimator::CheckSplineMeta(const SplineMeta<SplineOrder>& spline_meta)
	{
		return true;

		const int numKnots = trajectory_->numKnots();
		if (numKnots <= SplineOrder)
		{
			return false;
		}

		for (auto const& seg : spline_meta.segments)
		{
			size_t start_idx = trajectory_->computeTIndex(seg.t0 + 1e-9).second; // 开始节点的标号
			size_t end_idx = start_idx + SplineOrder;
			if (end_idx >= numKnots)
				return false;
		}

		return true;
	}

	// 插入控制节点
	void TrajectoryEstimator::AddControlPoints(const SplineMeta<SplineOrder>& spline_meta, std::vector<double*>& vec, bool addPosKont)
	{
		for (auto const& seg : spline_meta.segments)
		{
			size_t start_idx = trajectory_->computeTIndex(seg.t0 + 1e-9).second; // 开始节点的标号
			for (size_t i = start_idx; i < (start_idx + seg.NumParameters()); ++i) {
				//printf("i = %zd numKnots = %zd\n", i, trajectory_->numKnots());
				if (addPosKont)
				{
					vec.emplace_back(trajectory_->getKnotPos(i).data());                // 位置
					problem_->AddParameterBlock(trajectory_->getKnotPos(i).data(), 3);
					if (options_.lock_trajP)
						problem_->SetParameterBlockConstant(vec.back());
				}
				else
				{
					vec.emplace_back(trajectory_->getKnotSO3(i).data());                // 姿态
					problem_->AddParameterBlock(trajectory_->getKnotSO3(i).data(), 4, local_parameterization);
					if (options_.lock_trajR)
						problem_->SetParameterBlockConstant(vec.back());

					Eigen::Quaterniond qq = trajectory_->getKnotSO3(i).unit_quaternion();
					//printf("%10zd %10.3lf %9.6lf %9.6lf %9.6lf %9.6lf\n", i, trajectory_->getTime(i), qq.x(), qq.y(), qq.z(), qq.w());
				}

				//   if (options_.lock_traj)
				   //{
				//     problem_->SetParameterBlockConstant(vec.back());
				//   }
			}
		}
	}

	void TrajectoryEstimator::AddLoamMeasurement(const PointCorrespondence& pc,
		double weight, double huber_loss) {

		// step 1: 
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {pc.t_map, pc.t_map}, {pc.t_point, pc.t_point} }, spline_meta);

		// step 2: 
		using Functor = PointFeatureFactor;
		Functor* functor = new Functor(pc, spline_meta, weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

		// step 3: 
		cost_function->SetNumResiduals(1);
		for (int i = 0; i < spline_meta.NumParameters(); i++)	cost_function->AddParameterBlock(4); /// add so3 knots
		for (int i = 0; i < spline_meta.NumParameters(); i++)	cost_function->AddParameterBlock(3); /// add vec3 knots
		cost_function->AddParameterBlock(4);   // R_LtoI
		cost_function->AddParameterBlock(3);   // p_LinI

		// step 4: 
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);
		AddControlPoints(spline_meta, vec, true);
		vec.emplace_back(calib_param_->so3_LtoI.data());
		problem_->AddParameterBlock(calib_param_->so3_LtoI.data(), 4, local_parameterization);
		vec.emplace_back(calib_param_->p_LinI.data());

		// step 5:
		ceres::HuberLoss loss_function_(huber_loss);
		problem_->AddResidualBlock(cost_function, NULL, vec);
	}

	void TrajectoryEstimator::AddLiDARPoseMeasurement(const PoseData& pose_data,
		double rot_weight,
		double pos_weight) {
		// step 1: 
		SplineMeta<SplineOrder> spline_meta;
		if (options_.lock_t_offset) {
			double t_pose = pose_data.timestamp + trajectory_->GetTrajParam(pose_data.timestamp)->time_offset;
			if (t_pose < trajectory_->minTime() || t_pose >= trajectory_->maxTime())
				return;
			trajectory_->CaculateSplineMeta({ {t_pose, t_pose} }, spline_meta);
		}
		else {
			double t_min = pose_data.timestamp - options_.t_offset_padding;
			double t_max = pose_data.timestamp + options_.t_offset_padding;
			if (t_min < trajectory_->minTime() || t_max >= trajectory_->maxTime())
				return;
			trajectory_->CaculateSplineMeta({ {t_min, t_max} }, spline_meta);
		}

		// step 2: 
		using Functor = LiDARPoseFactor;
		Functor* functor = new Functor(pose_data, spline_meta, rot_weight, pos_weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

		// step 3: 
		cost_function->SetNumResiduals(6);
		for (int i = 0; i < spline_meta.NumParameters(); i++)	cost_function->AddParameterBlock(4); /// add so3 knots
		for (int i = 0; i < spline_meta.NumParameters(); i++)	cost_function->AddParameterBlock(3); /// add vec3 knots
		cost_function->AddParameterBlock(4);  // R_LtoI
		cost_function->AddParameterBlock(3);  // p_LinI
		cost_function->AddParameterBlock(1);  // time_offset

		// step 4: 
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);
		AddControlPoints(spline_meta, vec, true);
		vec.emplace_back(calib_param_->so3_LtoI.data());
		problem_->AddParameterBlock(calib_param_->so3_LtoI.data(), 4, local_parameterization);
		vec.emplace_back(calib_param_->p_LinI.data());

		auto param = trajectory_->GetTrajParam(pose_data.timestamp);
		vec.emplace_back(&param->time_offset);
		problem_->AddParameterBlock(&param->time_offset, 1);

		if (options_.lock_t_offset) {
			problem_->SetParameterBlockConstant(&param->time_offset);
		}
		else {
			problem_->SetParameterLowerBound(&param->time_offset, 0, -options_.t_offset_padding);
			problem_->SetParameterUpperBound(&param->time_offset, 0, options_.t_offset_padding);
		}

		// step 4: 
		problem_->AddResidualBlock(cost_function, new ceres::CauchyLoss(0.5), vec);
	}

	void TrajectoryEstimator::AddLiDARSurfelMeasurement(const PointCorrespondence& pc, double lidar_weight) {

		// step 1: 
		SplineMeta<SplineOrder> spline_meta;
		if (options_.lock_t_offset)
		{
			double t_lidar = pc.t_point + trajectory_->GetTrajParam(pc.t_point)->time_offset;
			if (t_lidar < trajectory_->minTime() || t_lidar >= trajectory_->maxTime())
				return;

			trajectory_->CaculateSplineMeta({ {t_lidar, t_lidar} }, spline_meta);
		}
		else
		{
			double t_min = pc.t_point - options_.t_offset_padding;
			double t_max = pc.t_point + options_.t_offset_padding;
			if (t_min < trajectory_->minTime() || t_max >= trajectory_->maxTime())
				return;
			trajectory_->CaculateSplineMeta({ {t_min, t_max} }, spline_meta);
		}

		// step 2: 
		using Functor = PointPlaneFactor;
		Functor* functor = new Functor(pc, spline_meta, lidar_weight, !options_.lock_LiDAR_intrinsic);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

		// step 3: 
		cost_function->SetNumResiduals(1);
		for (int i = 0; i < spline_meta.NumParameters(); i++)
			cost_function->AddParameterBlock(4); /// add so3 knots
		for (int i = 0; i < spline_meta.NumParameters(); i++)
			cost_function->AddParameterBlock(3); /// add vec3 knots
		cost_function->AddParameterBlock(4);  // R_LtoI
		cost_function->AddParameterBlock(3);  // p_LinI
		cost_function->AddParameterBlock(1);  // time_offset

		// step 4: 
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);
		AddControlPoints(spline_meta, vec, true);
		vec.emplace_back(calib_param_->so3_LtoI.data());
		problem_->AddParameterBlock(calib_param_->so3_LtoI.data(), 4, local_parameterization);
		vec.emplace_back(calib_param_->p_LinI.data());
		auto param = trajectory_->GetTrajParam(pc.t_point);
		vec.emplace_back(&param->time_offset);

		problem_->AddParameterBlock(&param->time_offset, 1);
		if (options_.lock_t_offset) {
			problem_->SetParameterBlockConstant(&param->time_offset);
		}
		else {
			problem_->SetParameterLowerBound(&param->time_offset, 0, -options_.t_offset_padding);
			problem_->SetParameterUpperBound(&param->time_offset, 0, options_.t_offset_padding);
		}

		if (!options_.lock_LiDAR_intrinsic) {
			/// opt_laser_param
			int laser_id = int(pc.point_raw.x());
			for (int idx = 0; idx < 6; ++idx) {
				cost_function->AddParameterBlock(1);  // laser_param
				auto laser_param = calib_param_->lidar_intrinsic.GetLaserParam(laser_id, idx);
				vec.emplace_back(laser_param);
				problem_->AddParameterBlock(laser_param, 1);

				// fix fist laser_id' param partly
				if (0 == laser_id && idx > 1) {
					problem_->SetParameterBlockConstant(laser_param);
				}
				// horiz_offset_mm
				//      if (idx == 3) {
				//        problem_->SetParameterBlockConstant(laser_param);
				//      }
				if (5 == idx) {
					problem_->SetParameterLowerBound(laser_param, 0, -0.017);  // 1 degree
					problem_->SetParameterUpperBound(laser_param, 0, 0.017);
				}
			}
		}

		// step 5:
		ceres::LossFunction* loss_function;
		loss_function = new ceres::CauchyLoss(1.0);  // adopt from vins-mono
		// problem_->AddResidualBlock(cost_function, NULL, vec);
		problem_->AddResidualBlock(cost_function, loss_function, vec);
	}


	/// @brief 添加陀螺数据
	///
	/// @param[in] imu_data	    IMU数据
	/// @param[in] gyro_weight	陀螺权重
	/// @return void
	void TrajectoryEstimator::AddIMUGyroMeasurement(const IMUData& imu_data, double gyro_weight)
	{
		// step 1: 计算当前时间的所有支撑节点
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {imu_data.timestamp, imu_data.timestamp} }, spline_meta);

		// step 2: 构建残差快
		using Functor = GyroscopeWithConstantBiasFactor;
		Functor* functor = new Functor(imu_data, spline_meta, gyro_weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

		// step 3: 设置残差块大小和参数块大小
		cost_function->SetNumResiduals(3);
		for (int i = 0; i < spline_meta.NumParameters(); i++)	cost_function->AddParameterBlock(4); /// add so3 knots
		cost_function->AddParameterBlock(3);  // gyro bias

		// step 4: 添加参数块
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);       // 姿态
		auto param = trajectory_->GetTrajParam(imu_data.timestamp); // 陀螺零偏
		vec.emplace_back(param->gyro_bias.data());
		if (options_.lock_wb)	problem_->SetParameterBlockConstant(param->gyro_bias.data());

		// step 5: 添加残差
		problem_->AddResidualBlock(cost_function, NULL, vec);

		return;
	}


	/// @brief 添加陀螺数据
	///
	/// @param[in] timestamp	当前时间
	/// @param[in] gyro_data	陀螺数据
	/// @param[in] gyro_weight	陀螺权重
	/// @return void
	void TrajectoryEstimator::AddIMUGyroMeasurement(const double timestamp, const Eigen::Vector3d& gyro_data, double gyro_weight)
	{
		// step 1
		IMUData imu_data;
		imu_data.timestamp = timestamp;
		imu_data.gyro = gyro_data;

		// step 2
		AddIMUGyroMeasurement(imu_data, gyro_weight);

		return;
	}


	/// @brief 添加加计数据
	///
	/// @param[in] imu_data	    IMU数据
	/// @param[in] accel_weight	加计权重
	/// @return void
	void TrajectoryEstimator::AddIMUAccelMeasurement(const IMUData& imu_data, double accel_weight)
	{
		// step 1: 计算当前时间的所有支撑节点
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {imu_data.timestamp, imu_data.timestamp} }, spline_meta);

		// 更新重力加速度
		options_.local_traj = false;
		Eigen::Vector3d gravity = options_.gravity;
		if (!options_.local_traj)
		{
			Eigen::Vector3d XYZ = trajectory_->position(imu_data.timestamp);
			Eigen::Vector3d LLH = XYZ2LLH(XYZ);
			gravity = getGravityECEF(LLH);
		}

		// step 2: 构建残差快
		using Functor = AccelerometerWithConstantBiasFactor;
		Functor* functor = new Functor(imu_data, spline_meta, accel_weight, options_.local_traj, gravity);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

		// step 3: 设置残差块大小和参数块大小
		cost_function->SetNumResiduals(3); // 残差3维
		for (int i = 0; i < spline_meta.NumParameters(); i++)	cost_function->AddParameterBlock(4); // 姿态
		for (int i = 0; i < spline_meta.NumParameters(); i++)	cost_function->AddParameterBlock(3); // 位置
		cost_function->AddParameterBlock(3);  // 加计零偏

		// step 4: 添加参数块
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);       // 姿态
		AddControlPoints(spline_meta, vec, true); // 位置
		auto param = trajectory_->GetTrajParam(imu_data.timestamp); // 根据时间选取对应时间段的重力和零偏
		vec.emplace_back(param->acce_bias.data()); problem_->AddParameterBlock(param->acce_bias.data(), 3); // 加计零偏
		if (options_.lock_ab)	problem_->SetParameterBlockConstant(param->acce_bias.data());

		// step 5: 添加残差
		problem_->AddResidualBlock(cost_function, NULL, vec);

		return;
	}

	/// @brief 添加加计数据
	///
	/// @param[in] timestamp	当前时间
	/// @param[in] accel_data	加计数据
	/// @param[in] accel_weight	加计权重
	/// @return void
	void TrajectoryEstimator::AddIMUAccelMeasurement(const double timestamp, const Eigen::Vector3d& accel_data, double accel_weight)
	{
		// step 1
		IMUData imu_data;
		imu_data.timestamp = timestamp;
		imu_data.gyro = accel_data;

		// step 2
		AddIMUGyroMeasurement(imu_data, accel_weight);

		return;
	}


	void TrajectoryEstimator::AddIMUMeasurement(const IMUData& imu_data, double gyro_weight, double accel_weight)
	{
		// step 1: 计算当前时间的所有支撑节点
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {imu_data.timestamp, imu_data.timestamp} }, spline_meta);

		// 更新重力加速度
		Eigen::Vector3d gravity = options_.gravity;
		if (!options_.local_traj)
		{
			Eigen::Vector3d XYZ = trajectory_->position(imu_data.timestamp);
			Eigen::Vector3d LLH = XYZ2LLH(XYZ);
			gravity = getGravityECEF(LLH);
		}

		// step 2: 构建残差快
		using Functor = GyroAcceWithConstantBiasFactor;
		Functor* functor = new Functor(imu_data, spline_meta, gyro_weight, accel_weight, options_.local_traj, gravity);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

		// step 3: 设置残差块大小和参数块大小
		cost_function->SetNumResiduals(6); // 残差6维
		for (int i = 0; i < spline_meta.NumParameters(); i++)	cost_function->AddParameterBlock(4); // 姿态
		for (int i = 0; i < spline_meta.NumParameters(); i++)	cost_function->AddParameterBlock(3); // 位置
		cost_function->AddParameterBlock(3);  // 陀螺零偏
		cost_function->AddParameterBlock(3);  // 加计零偏
		cost_function->AddParameterBlock(2);  // 重力

		// step 4: 添加参数块
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);       // 姿态
		AddControlPoints(spline_meta, vec, true); // 位置
		auto param = trajectory_->GetTrajParam(imu_data.timestamp); // 根据时间选取对应时间段的重力和零偏
		vec.emplace_back(param->gyro_bias.data()); problem_->AddParameterBlock(param->gyro_bias.data(), 3); // 陀螺零偏
		if (options_.lock_wb)	problem_->SetParameterBlockConstant(param->gyro_bias.data());
		vec.emplace_back(param->acce_bias.data()); problem_->AddParameterBlock(param->acce_bias.data(), 3); // 加计零偏
		if (options_.lock_ab)	problem_->SetParameterBlockConstant(param->acce_bias.data());
		vec.emplace_back(param->g_refine.data());  problem_->AddParameterBlock(param->g_refine.data(), 2);  // 重力
		if (options_.lock_g)	problem_->SetParameterBlockConstant(param->g_refine.data());
		// for (int i = 0; i < 2; i++) {
		//   problem_->SetParameterLowerBound(param->g_refine.data(), i, -M_PI / 2.0);
		//   problem_->SetParameterUpperBound(param->g_refine.data(), i, M_PI / 2.0);
		// }

		/// opt_imu_param
		{
			cost_function->AddParameterBlock(6);  // 陀螺尺度因子(0-2)和非正交误差(3-5)
			cost_function->AddParameterBlock(6);  // 加计尺度因子(0-2)和非正交误差(3-5)
			cost_function->AddParameterBlock(9);  // 加计对陀螺的影响
			cost_function->AddParameterBlock(4);  // 陀螺相对于加计的安置角

			auto M_w_param = calib_param_->imu_intrinsic.GetMwVector();
			auto M_a_param = calib_param_->imu_intrinsic.GetMaVector();
			auto A_w_param = calib_param_->imu_intrinsic.GetAwVector();
			auto q_WtoA_param = calib_param_->imu_intrinsic.GetQWtoAVector();
			vec.emplace_back(M_w_param);
			vec.emplace_back(M_a_param);
			vec.emplace_back(A_w_param);
			vec.emplace_back(q_WtoA_param);

			problem_->AddParameterBlock(M_w_param, 6);
			problem_->AddParameterBlock(M_a_param, 6);
			problem_->AddParameterBlock(A_w_param, 9);
			problem_->AddParameterBlock(q_WtoA_param, 4, local_parameterization);

			if (options_.lock_IMU_intrinsic) {
				problem_->SetParameterBlockConstant(M_w_param);
				problem_->SetParameterBlockConstant(M_a_param);
				problem_->SetParameterBlockConstant(A_w_param);
				problem_->SetParameterBlockConstant(q_WtoA_param);
			}
			else {
				// 不考虑 G sensitivity
				problem_->SetParameterBlockConstant(A_w_param);
			}
			
			if (false)
			{
				// Subset parametrization fix partly param in SubsetParameterization
				/// 0 3 6
				/// 1 4 7
				/// 2 5 8
				static std::vector<int> vec_constant_param = { 1, 2, 5 };
				// vec_constant_param.insert(vec_constant_param.end(), {1,2,5});
				static ceres::SubsetParameterization* subset_parameterization = new ceres::SubsetParameterization(9, vec_constant_param);
				// problem_->AddParameterBlock(A_w_param, 9, subset_parameterization);
				problem_->AddParameterBlock(A_w_param, 9);
				problem_->SetParameterization(A_w_param, subset_parameterization);
			}
		}

		problem_->AddResidualBlock(cost_function, NULL, vec);
	}


	void TrajectoryEstimator::AddPoseMeasurement(const PoseData& pose_data, double rot_weight, double pos_weight) {
		// step 1: 计算节点参数
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {pose_data.timestamp, pose_data.timestamp} }, spline_meta);

		// step 2: 
		using Functor = IMUPoseFactor;
		Functor* functor = new Functor(pose_data, spline_meta, rot_weight, pos_weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

		// step 3: 
		cost_function->SetNumResiduals(6);
		for (int i = 0; i < spline_meta.NumParameters(); i++)
			cost_function->AddParameterBlock(4); /// add so3 knots
		for (int i = 0; i < spline_meta.NumParameters(); i++)
			cost_function->AddParameterBlock(3); /// add vec3 knots

		  // step 4: 
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);
		AddControlPoints(spline_meta, vec, true);
		problem_->AddResidualBlock(cost_function, NULL, vec);
	}

	/// @brief IMU位置
	///
	/// @param[in] pose_data	IMU位姿
	/// @param[in] pos_weight	位置权重
	/// @return void
	void TrajectoryEstimator::AddPositionMeasurement(const PoseData& pose_data, double pos_weight)
	{
		// step 1: 获取当前时间的支撑节点
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {pose_data.timestamp, pose_data.timestamp} }, spline_meta);

		// step 2: 构建残差块
		using Functor = IMUPositionFactor;
		Functor* functor = new Functor(pose_data, spline_meta, pos_weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

		// step 3: 设置残差块大小和参数块大小
		cost_function->SetNumResiduals(3);
		for (int i = 0; i < spline_meta.NumParameters(); i++)
			cost_function->AddParameterBlock(3);

		// step 4: 添加参数块
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec, true);
		problem_->AddResidualBlock(cost_function, NULL, vec);
	}

	/// @brief IMU位置
	///
	/// @param[in] timestamp	当前时间
	/// @param[in] position	IMU位置
	/// @param[in] pos_weight	位置权重
	/// @return void
	void TrajectoryEstimator::AddPositionMeasurement(const double timestamp, const Eigen::Vector3d position, double pos_weight)
	{
		// step 1: 获取当前时间的支撑节点
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {timestamp, timestamp} }, spline_meta);

		// step 2: 构建残差块
		using Functor = IMUPositionFactor;
		Functor* functor = new Functor(timestamp, position, spline_meta, pos_weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

		// step 3: 设置残差块大小和参数块大小
		cost_function->SetNumResiduals(3);
		for (int i = 0; i < spline_meta.NumParameters(); i++)	cost_function->AddParameterBlock(3);

		// step 4: 添加参数块
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec, true);

		// step 5: 添加残差
		problem_->AddResidualBlock(cost_function, NULL, vec);
	}

	/// @brief IMU速度
	///
	/// @param[in] pose_data	IMU位姿
	/// @param[in] pos_weight	速度权重
	/// @return void
	void TrajectoryEstimator::AddVelocityMeasurement(const PoseData& pose_data, double vel_weight)
	{
		// step 1: 获取当前时间的支撑节点
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {pose_data.timestamp, pose_data.timestamp} }, spline_meta);

		// step 2: 构建残差块
		using Functor = IMUVelocityFactor;
		Functor* functor = new Functor(pose_data, spline_meta, vel_weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

		// step 3: 设置残差块大小和参数块大小
		cost_function->SetNumResiduals(3);
		for (int i = 0; i < spline_meta.NumParameters(); i++)	cost_function->AddParameterBlock(3);

		// step 4: 添加参数块
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec, true);

		// step 5: 添加残差
		problem_->AddResidualBlock(cost_function, NULL, vec);
	}

	/// @brief IMU速度
	///
	/// @param[in] timestamp	当前时间
	/// @param[in] position	IMU速度
	/// @param[in] pos_weight	速度权重
	/// @return void
	void TrajectoryEstimator::AddVelocityMeasurement(const double timestamp, const Eigen::Vector3d Velocity, double vel_weight)
	{
		// step 1: 获取当前时间的支撑节点
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {timestamp, timestamp} }, spline_meta);

		// step 2: 构建残差块
		using Functor = IMUVelocityFactor;
		Functor* functor = new Functor(timestamp, Velocity, spline_meta, vel_weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

		// step 3: 设置残差块大小和参数块大小
		cost_function->SetNumResiduals(3);
		for (int i = 0; i < spline_meta.NumParameters(); i++)
			cost_function->AddParameterBlock(3);

		// step 4: 添加参数块
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec, true);

		// step 5: 添加残差
		problem_->AddResidualBlock(cost_function, NULL, vec);
	}

	/// @brief IMU姿态
	///
	/// @param[in] pose_data	IMU位姿
	/// @param[in] pos_weight	姿态权重
	/// @return void
	void TrajectoryEstimator::AddOrientationMeasurement(const PoseData& pose_data, double rot_weight)
	{
		// step 1: 获取当前时间的支撑节点
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {pose_data.timestamp, pose_data.timestamp} }, spline_meta);

		// step 2: 构建残差块
		using Functor = IMUOrientationFactor;
		Functor* functor = new Functor(pose_data, spline_meta, rot_weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

		// step 3: 设置残差块大小和参数块大小
		cost_function->SetNumResiduals(3);
		for (int i = 0; i < spline_meta.NumParameters(); i++)	cost_function->AddParameterBlock(4);

		// step 4: 添加参数块
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);

		// step 5: 添加残差
		problem_->AddResidualBlock(cost_function, NULL, vec);
	}

	/// @brief IMU姿态
	///
	/// @param[in] timestamp   当前时间
	/// @param[in] Orientation IMU姿态
	/// @param[in] rot_weight  姿态权重
	/// @return void
	void TrajectoryEstimator::AddOrientationMeasurement(const double timestamp, const SO3d& Orientation, double rot_weight)
	{
		// step 1: 获取当前时间的支撑节点
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {timestamp, timestamp} }, spline_meta);

		// step 2: 构建残差块
		using Functor = IMUOrientationFactor;
		Functor* functor = new Functor(timestamp, Orientation, spline_meta, rot_weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

		// step 3: 设置残差块大小和参数块大小
		cost_function->SetNumResiduals(3);
		for (int i = 0; i < spline_meta.NumParameters(); i++)	cost_function->AddParameterBlock(4);

		// step 4: 添加参数块
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);

		// step 5: 添加残差
		problem_->AddResidualBlock(cost_function, NULL, vec);

	}

	void TrajectoryEstimator::AddAccelBiasPriorConstraint(const Eigen::Vector3d weight)
	{
		if (calib_param_->segment_param.size() < 1)
			return;

		for (size_t i = 0; i < calib_param_->segment_param.size(); i++)
		{
			// step 1 : 
			auto param = calib_param_->segment_param[i].acce_bias.data();
			if (!problem_->HasParameterBlock(param))
				continue;

			// step 2: 设置残差块大小和参数块大小
			using Functor = IMUBiasPriorFactor;
			Functor* functor = new Functor(weight);
			auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

			// step 3: 设置残差块大小和参数块大小
			cost_function->SetNumResiduals(3);
			cost_function->AddParameterBlock(3);

			// step 4: 添加参数块
			std::vector<double*> vec;
			vec.emplace_back(calib_param_->segment_param[i].acce_bias.data());
			problem_->AddResidualBlock(cost_function, NULL, vec);
		}

		return;
	}
	void TrajectoryEstimator::AddAccelBiasPriorConstraint(const double weight)
	{
		Eigen::Vector3d	weight_;
		weight_(0) = weight_(1) = weight_(2) = weight;
		return AddAccelBiasPriorConstraint(weight_);
	}

	void TrajectoryEstimator::AddAccelBiasPriorConstraint(const Eigen::Vector3d &bias, const Eigen::Vector3d weight)
	{
		if (calib_param_->segment_param.size() < 1)
			return;

		for (size_t i = 0; i < calib_param_->segment_param.size(); i++)
		{
			// step 1 : 
			auto param = calib_param_->segment_param[i].acce_bias.data();
			if (!problem_->HasParameterBlock(param))
				continue;

			// step 2: 设置残差块大小和参数块大小
			using Functor = IMUBiasPriorFactor;
			Functor* functor = new Functor(bias, weight);
			auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

			// step 3: 设置残差块大小和参数块大小
			cost_function->SetNumResiduals(3);
			cost_function->AddParameterBlock(3);

			// step 4: 添加参数块
			std::vector<double*> vec;
			vec.emplace_back(param);
			problem_->AddResidualBlock(cost_function, NULL, vec);
		}

		return;
	}
	void TrajectoryEstimator::AddAccelBiasPriorConstraint(const Eigen::Vector3d &bias, const double weight)
	{
		Eigen::Vector3d	weight_;
		weight_(0) = weight_(1) = weight_(2) = weight;
		return AddAccelBiasPriorConstraint(bias, weight_);
	}
	void TrajectoryEstimator::AddAccelBiasPriorConstraint(const double bias, const double weight)
	{
		Eigen::Vector3d	bias_;
		Eigen::Vector3d	weight_;
		bias_(0) = bias_(1) = bias_(2) = bias;
		weight_(0) = weight_(1) = weight_(2) = weight;
		return AddAccelBiasPriorConstraint(bias_, weight_);
	}

	void TrajectoryEstimator::AddGyroBiasPriorConstraint(const Eigen::Vector3d weight)
	{
		if (calib_param_->segment_param.size() < 1)
			return;

		for (size_t i = 0; i < calib_param_->segment_param.size(); i++)
		{
			// step 1 : 
			auto param = calib_param_->segment_param[i].gyro_bias.data();
			if (!problem_->HasParameterBlock(param))
				continue;

			// step 2: 设置残差块大小和参数块大小
			using Functor = IMUBiasPriorFactor;
			Functor* functor = new Functor(weight);
			auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

			// step 3: 设置残差块大小和参数块大小
			cost_function->SetNumResiduals(3);
			cost_function->AddParameterBlock(3);

			// step 4: 添加参数块
			std::vector<double*> vec;
			vec.emplace_back(param);
			problem_->AddResidualBlock(cost_function, NULL, vec);
		}

		return;
	}
	void TrajectoryEstimator::AddGyroBiasPriorConstraint(const double weight)
	{
		Eigen::Vector3d	weight_;
		weight_(0) = weight_(1) = weight_(2) = weight;
		return AddGyroBiasPriorConstraint(weight_);
	}
	void TrajectoryEstimator::AddGyroBiasPriorConstraint(const Eigen::Vector3d &bias, const Eigen::Vector3d weight)
	{
		if (calib_param_->segment_param.size() < 1)
			return;

		for (size_t i = 0; i < calib_param_->segment_param.size(); i++)
		{
			// step 1 : 
			auto param = calib_param_->segment_param[i].gyro_bias.data();
			if (!problem_->HasParameterBlock(param))
				continue;

			// step 2: 设置残差块大小和参数块大小
			using Functor = IMUBiasPriorFactor;
			Functor* functor = new Functor(bias, weight);
			auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

			// step 3: 设置残差块大小和参数块大小
			cost_function->SetNumResiduals(3);
			cost_function->AddParameterBlock(3);

			// step 4: 添加参数块
			std::vector<double*> vec;
			vec.emplace_back(param);
			problem_->AddResidualBlock(cost_function, NULL, vec);
		}

		return;
	}
	void TrajectoryEstimator::AddGyroBiasPriorConstraint(const Eigen::Vector3d &bias, const double weight)
	{
		Eigen::Vector3d	weight_;
		weight_(0) = weight_(1) = weight_(2) = weight;
		return AddGyroBiasPriorConstraint(bias, weight_);
	}
	void TrajectoryEstimator::AddGyroBiasPriorConstraint(const double bias, const double weight)
	{
		Eigen::Vector3d	bias_;
		Eigen::Vector3d	weight_;
		bias_(0) = bias_(1) = bias_(2) = bias;
		weight_(0) = weight_(1) = weight_(2) = weight;
		return AddGyroBiasPriorConstraint(bias_, weight_);
	}

	void TrajectoryEstimator::AddAccelPriorBias(const double bias, const double biasStd)
	{
		if (!calib_param_ || calib_param_->segment_param.size() < 1)
			return;

		double	weight = (1.0 / (biasStd));
		for (size_t i = 0; i < calib_param_->segment_param.size(); i++)
		{
			// step 0 : 检查参数块
			auto param = calib_param_->segment_param[i].acce_bias.data();
			if (!problem_->HasParameterBlock(param))
				return;

			// step 1: 
			using Functor = IMUBiasPriorFactor;
			Functor* functor = new Functor(bias, weight);
			auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

			// step 2: 设置残差块大小和参数块大小
			cost_function->SetNumResiduals(3);
			cost_function->AddParameterBlock(3); // 前一个

			// step 3: 添加参数块
			std::vector<double*> vec;
			vec.emplace_back(param);

			// step 4: 添加残差
			problem_->AddResidualBlock(cost_function, NULL, vec);
		}

		return;
	}
	void TrajectoryEstimator::AddAccelPriorBias(const double biasStd)
	{
		const double bias = 0.0;
		AddAccelPriorBias(bias, biasStd);
		return;
	}

	void TrajectoryEstimator::AddGyroPriorBias(const double bias, const double biasStd)
	{
		if (!calib_param_ || calib_param_->segment_param.size() < 1)
			return;

		double	weight = sqrt(1.0 / (biasStd));
		for (size_t i = 0; i < calib_param_->segment_param.size(); i++)
		{
			// step 0 : 检查参数块
			auto param = calib_param_->segment_param[i].gyro_bias.data();
			if (!problem_->HasParameterBlock(param))
				return;

			// step 1: 
			using Functor = IMUBiasPriorFactor;
			Functor* functor = new Functor(bias, weight);
			auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

			// step 2: 设置残差块大小和参数块大小
			cost_function->SetNumResiduals(3);
			cost_function->AddParameterBlock(3); // 前一个

			// step 3: 添加参数块
			std::vector<double*> vec;
			vec.emplace_back(param);

			// step 4: 添加残差
			problem_->AddResidualBlock(cost_function, NULL, vec);
		}

		return;

	}
	void TrajectoryEstimator::AddGyroPriorBias(const double biasStd)
	{
		const double bias = 0.0;
		AddGyroPriorBias(bias, biasStd);
		return;
	}

	/// @brief 加计零偏随机游走约束
	///
	/// @param[in] PNSD				加计零偏随机游走谱密度
	/// @return void
	void TrajectoryEstimator::AddAccelBiasBetweenConstraint(const Eigen::Vector3d PNSD)
	{
		if (calib_param_->segment_param.size() < 2)
			return;

		// 相邻弧段之间的约束
		for (size_t i = 1; i < calib_param_->segment_param.size(); i++)
		{
			// step 0 : 检查参数块
			auto param1 = calib_param_->segment_param[i - 1].acce_bias.data();
			auto param2 = calib_param_->segment_param[i].acce_bias.data();
			if (!problem_->HasParameterBlock(param1) || !problem_->HasParameterBlock(param2))
				continue;
			
			if (calib_param_->segment_param[i - 1].starttime == 0.0 || calib_param_->segment_param[i - 1].endtime == 0.0 ||
				calib_param_->segment_param[i].starttime == 0.0 || calib_param_->segment_param[i].endtime == 0.0)
				continue;
			if (calib_param_->segment_param[i - 1].starttime >= calib_param_->segment_param[i - 1].endtime ||
				calib_param_->segment_param[i].starttime >= calib_param_->segment_param[i].endtime)
				continue;

			// 取中间时间
			double	t1 = (calib_param_->segment_param[i - 1].starttime + calib_param_->segment_param[i - 1].endtime) / 2;
			double	t2 = (calib_param_->segment_param[i].starttime + calib_param_->segment_param[i].endtime) / 2;
			double	dt = t2 - t1;
			Eigen::Vector3d	weight_;
			weight_(0) = 1.0 / (PNSD(0) * sqrt(dt));
			weight_(1) = 1.0 / (PNSD(1) * sqrt(dt));
			weight_(2) = 1.0 / (PNSD(2) * sqrt(dt));
			//weight_(0) = sqrt(1.0 / (SQR(PNSD(0)) * dt));
			//weight_(1) = sqrt(1.0 / (SQR(PNSD(1)) * dt));
			//weight_(2) = sqrt(1.0 / (SQR(PNSD(2)) * dt));
			//weight_(0) = 1.0 / (PNSD(0) * sqrt(dt));
			//weight_(1) = 1.0 / (PNSD(1) * sqrt(dt));
			//weight_(2) = 1.0 / (PNSD(2) * sqrt(dt));

			// step 1: 
			using Functor = IMUBiasBetweenFactor;
			Functor* functor = new Functor(weight_);
			auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

			// step 2: 设置残差块大小和参数块大小
			cost_function->SetNumResiduals(3);
			cost_function->AddParameterBlock(3); // 前一个
			cost_function->AddParameterBlock(3); // 后一个

			// step 3: 添加参数块
			std::vector<double*> vec;
			vec.emplace_back(param1); vec.emplace_back(param2);

			// step 4: 添加残差
			problem_->AddResidualBlock(cost_function, NULL, vec);
		}

		return;
	}

	/// @brief 加计零偏随机游走约束
	///
	/// @param[in] PNSD				加计零偏随机游走谱密度
	/// @return void
	void TrajectoryEstimator::AddAccelBiasBetweenConstraint(const double PNSD)
	{
		Eigen::Vector3d PNSD_;
		PNSD_(0) = PNSD_(1) = PNSD_(2) = PNSD;
		return AddAccelBiasBetweenConstraint(PNSD_);
	}

	/// @brief 陀螺零偏随机游走约束
	///
	/// @param[in] PNSD				陀螺零偏随机游走谱密度
	/// @return void
	void TrajectoryEstimator::AddGyroBiasBetweenConstraint(const Eigen::Vector3d PNSD)
	{
		if (calib_param_->segment_param.size() < 2)
			return;

		// 相邻弧段之间的约束
		for (size_t i = 1; i < calib_param_->segment_param.size(); i++)
		{
			// step 1 : 
			auto param1 = calib_param_->segment_param[i - 1].gyro_bias.data();
			auto param2 = calib_param_->segment_param[i].gyro_bias.data();
			if (!problem_->HasParameterBlock(param1) || !problem_->HasParameterBlock(param2))
				continue;

			if (calib_param_->segment_param[i - 1].starttime == 0.0 || calib_param_->segment_param[i - 1].endtime == 0.0 ||
				calib_param_->segment_param[i].starttime == 0.0 || calib_param_->segment_param[i].endtime == 0.0)
				continue;
			if (calib_param_->segment_param[i - 1].starttime >= calib_param_->segment_param[i - 1].endtime ||
				calib_param_->segment_param[i].starttime >= calib_param_->segment_param[i].endtime)
				continue;

			double	t1 = (calib_param_->segment_param[i - 1].starttime + calib_param_->segment_param[i - 1].endtime) / 2;
			double	t2 = (calib_param_->segment_param[i].starttime + calib_param_->segment_param[i].endtime) / 2;
			double	dt = t2 - t1;
			Eigen::Vector3d	weight_;
			weight_(0) = sqrt(1.0 / (PNSD(0) * dt));
			weight_(1) = sqrt(1.0 / (PNSD(1) * dt));
			weight_(2) = sqrt(1.0 / (PNSD(2) * dt));

			// step 2: 设置残差块大小和参数块大小
			using Functor = IMUBiasBetweenFactor;
			Functor* functor = new Functor(weight_);
			auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

			// step 3: 设置残差块大小和参数块大小
			cost_function->SetNumResiduals(3);
			cost_function->AddParameterBlock(3); cost_function->AddParameterBlock(3);

			// step 4: 添加参数块
			std::vector<double*> vec;
			vec.emplace_back(param1); vec.emplace_back(param2);
			problem_->AddResidualBlock(cost_function, NULL, vec);
		}

		return;
	}

	/// @brief 陀螺零偏随机游走约束
	///
	/// @param[in] PNSD				陀螺零偏随机游走谱密度
	/// @return void
	void TrajectoryEstimator::AddGyroBiasBetweenConstraint(const double PNSD)
	{
		Eigen::Vector3d PNSD_;
		PNSD_(0) = PNSD_(1) = PNSD_(2) = PNSD;
		return AddGyroBiasBetweenConstraint(PNSD_);
	}

	/// @brief 零偏随机游走约束
	///
	/// @param[in] GyrPNSD			陀螺零偏随机游走谱密度
	/// @param[in] AccPNSD			加计零偏随机游走谱密度
	/// @return void
	void TrajectoryEstimator::AddIMUBiasBetweenConstraint(const Eigen::Vector3d GyrPNSD, const Eigen::Vector3d AccPNSD)
	{
		AddGyroBiasBetweenConstraint(GyrPNSD);
		AddAccelBiasBetweenConstraint(AccPNSD);

	}

	/// @brief 零偏随机游走约束
	///
	/// @param[in] GyrPNSD			陀螺零偏随机游走谱密度
	/// @param[in] AccPNSD			加计零偏随机游走谱密度
	/// @return void
	void TrajectoryEstimator::AddIMUBiasBetweenConstraint(const double GyrPNSD, const double AccPNSD)
	{
		AddGyroBiasBetweenConstraint(GyrPNSD);
		AddAccelBiasBetweenConstraint(AccPNSD);
	}


	void TrajectoryEstimator::AddQuadraticIntegralFactor(double min_time,
		double max_time,
		Eigen::Matrix3d weight) {
		size_t min_s = trajectory_->computeTIndex(min_time).second;
		size_t max_s = trajectory_->computeTIndex(max_time).second;

		std::vector<double> times;
		for (size_t i = min_s; i <= max_s; i++) {
			double timestamp = trajectory_->minTime() + i * trajectory_->getDt();
			times.push_back(timestamp);
		}
		times.back() += 1e-10;
		std::cout << YELLOW << "[AddQuadraticIntegralFactor] [" << std::fixed
			<< std::setprecision(5)
			<< trajectory_->minTime() + min_s * trajectory_->getDt() << ", "
			<< trajectory_->minTime() + max_s * trajectory_->getDt() << "]\n"
			<< RESET;

		if (times.empty()) return;

		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {times.front(), times.back()} }, spline_meta);

		using Functor = QuadraticIntegralFactor<SplineOrder, 2>;
		Functor* functor = new Functor(times, spline_meta, weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

		/// add R3 knots
		for (int i = 0; i < spline_meta.NumParameters(); i++) {
			cost_function->AddParameterBlock(3);
		}

		cost_function->SetNumResiduals(3 * SplineOrder * times.size());

		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec, true);

		problem_->AddResidualBlock(cost_function, NULL, vec);
	}

	void TrajectoryEstimator::AddAngularVelocityConvexHullFactor(double time, double weight) {
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {time, time} }, spline_meta);

		using Functor = AngularVelocityConvexHullFactor;
		Functor* functor = new Functor(spline_meta, weight);
		auto* cost_function =
			new ceres::DynamicAutoDiffCostFunction<Functor>(functor);

		/// add so3 knots
		for (int i = 0; i < spline_meta.NumParameters(); i++) {
			cost_function->AddParameterBlock(4);
		}

		size_t v_kont_num = spline_meta.NumParameters() - spline_meta.segments.size();
		cost_function->SetNumResiduals(v_kont_num);

		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);

		problem_->AddResidualBlock(cost_function, NULL, vec);
	}

	void TrajectoryEstimator::printPointResiduals()
	{
		if (pFrameContainer == NULL || pFrameContainer->MapPoints.size() <= 0)
			return;

		// 获取所有残差块
		std::vector<ceres::ResidualBlockId> residual_blocks;
		problem_->GetResidualBlocks(&residual_blocks);

		// 遍历所有地图点
		std::map<int64_t, MapPoint*>::iterator iter1;
		std::list<std::pair<ceres::ResidualBlockId, ObsPair>>::iterator iter2;
		for (iter1 = pFrameContainer->MapPoints.begin(); iter1 != pFrameContainer->MapPoints.end(); iter1++)
		{
			// 遍历所有残差块
			double	cost = 0.0, residuals[6] = { 0.0 };
			for (iter2 = iter1->second->residuals_blkid.begin(); iter2 != iter1->second->residuals_blkid.end(); iter2++)
			{
				ceres::ResidualBlockId	blkId = iter2->first;
				const int num_residuals = problem_->GetCostFunctionForResidualBlock(blkId)->num_residuals();
				problem_->EvaluateResidualBlock(blkId, false, &cost, residuals, NULL);

				// 根据残差个数自适应输出
				fprintf(stdout, "%p", blkId);
				for (int ires = 0; ires < num_residuals; ires++)
					fprintf(stdout, " %10.6lf", residuals[ires]);
				fprintf(stdout, " %10.6lf\n", cost);
			}
		}

		return;
	}

	bool TrajectoryEstimator::AddPointMeasurements(const LiDMeasType type, const double weight)
	{
		if (pFrameContainer == NULL || pFrameContainer->MapPoints.size() <= 0)
			return false;

		// 遍历所有地图点
		int		nPoint = 0, nObs = 0;
		std::map<int64_t, MapPoint*>::iterator iter;
		for (iter = pFrameContainer->MapPoints.begin(); iter != pFrameContainer->MapPoints.end(); ++iter)
		{
			iter->second->residuals_blkid.clear();
			if (iter->second->m_pFrames.size() < 2)
				continue;

			std::vector<std::pair<lpostk::Frame*, size_t>> vec(iter->second->m_pFrames.begin(), iter->second->m_pFrames.end());
			std::sort(vec.begin(), vec.end(), CompareFrameIdPair);
			nPoint++;

			int		nObs0 = 0;
			if (type == LiDMeasType::LOCALTOGLOBAL) // 每一帧与全局系
			{
				for (size_t i = 0; i < vec.size(); i++)
				{
					size_t innerId = vec[i].second;
					PointObs obs = vec[i].first->Points[innerId];
					obs.timestamp = fmod(obs.timestamp, 604800.0);
					ceres::ResidualBlockId blkId = AddPointLocalToGlobalMeasurement(obs, iter->second->XYZ, weight);

					ObsPair	obspair;
					obspair.Id1 = innerId; obspair.pFrame1 = vec[i].first;
					iter->second->residuals_blkid.push_back(std::pair<ceres::ResidualBlockId, ObsPair>(blkId, obspair));
					nObs0++;
				}
			}
			else if (type == LiDMeasType::LOCALTOLOCAL) // 任意两帧之间
			{
				for (size_t i = 0; i < vec.size() - 1; i++)
				{
					// 第一帧
					size_t innerId1 = vec[i].second;
					PointObs obs1 = vec[i].first->Points[innerId1];

					// 第二帧
					for (size_t j = i + 1; j < vec.size(); j++)
					{
						size_t innerId2 = vec[j].second;
						PointObs obs2 = vec[j].first->Points[innerId2];

						// 只保留周秒
						obs1.timestamp = fmod(obs1.timestamp, 604800.0);
						obs2.timestamp = fmod(obs2.timestamp, 604800.0);

						//if (j - i < 20)	continue;
						ceres::ResidualBlockId blkId = AddPointLocalToLocalMeasurement(obs1, obs2, weight);

						ObsPair	obspair;
						obspair.Id1 = innerId1; obspair.Id2 = innerId2;
						obspair.pFrame1 = vec[i].first; obspair.pFrame2 = vec[j].first;
						iter->second->residuals_blkid.push_back(std::pair<ceres::ResidualBlockId, ObsPair>(blkId, obspair));
						nObs0++;
					}
				}

				fprintf(stdout, "AddPointMeasurements ID = %010lld nFrame = %8zd nObs = %d\n", iter->first, iter->second->m_pFrames.size(), nObs0);
			}

			nObs += nObs0;
		}

		fprintf(stdout, "AddPointMeasurements nPoint = %d nObs = %d\n", nPoint, nObs);
		return nPoint > 0;
	}

	ceres::ResidualBlockId TrajectoryEstimator::AddPointLocalToGlobalMeasurement(const PointObs &obs, double PW[3], double weight)
	{
		// 计算样条节点位置
		const double time = obs.timestamp;
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {time, time} }, spline_meta);

		PointLocalToGlobalFactor* functor = new PointLocalToGlobalFactor(obs, spline_meta, weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<PointLocalToGlobalFactor>(functor);
		cost_function->SetNumResiduals(3);

		// 参数块大小
		for (int i = 0; i < spline_meta.NumParameters(); i++) 	cost_function->AddParameterBlock(4); // 姿态
		for (int i = 0; i < spline_meta.NumParameters(); i++)	cost_function->AddParameterBlock(3); // 位置
		cost_function->AddParameterBlock(4); // 外参R
		cost_function->AddParameterBlock(3); // 外参T
		cost_function->AddParameterBlock(3); // 陆标点

		// 参数块列表
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);       // 姿态
		AddControlPoints(spline_meta, vec, true); // 位置
		vec.emplace_back(calib_param_->so3_LtoI.data()); problem_->AddParameterBlock(calib_param_->so3_LtoI.data(), 4, local_parameterization); // 外参R
		vec.emplace_back(calib_param_->p_LinI.data()); problem_->AddParameterBlock(calib_param_->p_LinI.data(), 3); // 外参T
		vec.emplace_back(PW); problem_->AddParameterBlock(PW, 3); // 陆标点

		if (0)
		{
			double	*residuals = NULL, **jacobians = NULL;
			const int num_residuals = cost_function->num_residuals();
			const std::vector<int32_t> block_sizes = cost_function->parameter_block_sizes();

			// 分配内存
			residuals = (double *)calloc(num_residuals, sizeof(double));
			jacobians = (double **)calloc(block_sizes.size(), sizeof(double*)); // 第一维索引参数块个数
			for (size_t i = 0; i < block_sizes.size(); i++)
				jacobians[i] = (double *)calloc(block_sizes[i] * num_residuals, sizeof(double)); // 第二位一次为每个残差对参数块内每个参数的偏导

			// 计算残差和雅克比
			cost_function->Evaluate(vec.data(), residuals, jacobians);

			// 释放每次
			free(residuals); residuals = NULL;
			for (size_t i = 0; i < block_sizes.size(); i++)
			{
				free(jacobians[i]); jacobians[i] = NULL;
			}
			jacobians = NULL;
		}

		// 添加残差块
		ceres::ResidualBlockId blkId = problem_->AddResidualBlock(cost_function, NULL, vec);
		return blkId;
	}

	void TrajectoryEstimator::printLineResiduals()
	{

	}

	ceres::ResidualBlockId TrajectoryEstimator::AddPointLocalToLocalMeasurement(const PointObs &obs1, const PointObs &obs2, double weight)
	{
		// 计算样条节点位置
		const double time1 = obs1.timestamp;
		const double time2 = obs2.timestamp;
		double	t1 = time1 <= time2 ? time1 : time2;
		double	t2 = time1 >= time2 ? time1 : time2;
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {t1, t2} }, spline_meta);

		PointLocalToLocalFactor* functor = new PointLocalToLocalFactor(obs1, obs2, spline_meta, weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<PointLocalToLocalFactor>(functor);
		cost_function->SetNumResiduals(6);

		// 参数块大小
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 姿态
			cost_function->AddParameterBlock(4);
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 位置
			cost_function->AddParameterBlock(3);
		cost_function->AddParameterBlock(4); // 外参R
		cost_function->AddParameterBlock(3); // 外参T

		// 参数块列表
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);       // 第一帧姿态
		AddControlPoints(spline_meta, vec, true); // 第一帧位置
		vec.emplace_back(calib_param_->so3_LtoI.data()); problem_->AddParameterBlock(calib_param_->so3_LtoI.data(), 4, local_parameterization); // 外参R
		vec.emplace_back(calib_param_->p_LinI.data()); problem_->AddParameterBlock(calib_param_->p_LinI.data(), 3); // 外参T

		if (0)
		{
			double	*residuals = NULL, **jacobians = NULL;
			const int num_residuals = cost_function->num_residuals();
			const std::vector<int32_t> block_sizes = cost_function->parameter_block_sizes();

			// 分配内存
			residuals = (double *)calloc(num_residuals, sizeof(double));
			jacobians = (double **)calloc(block_sizes.size(), sizeof(double*)); // 第一维索引参数块个数
			for (size_t i = 0; i < block_sizes.size(); i++)
				jacobians[i] = (double *)calloc(block_sizes[i] * num_residuals, sizeof(double)); // 第二位一次为每个残差对参数块内每个参数的偏导

			// 计算残差和雅克比
			cost_function->Evaluate(vec.data(), residuals, jacobians);

			// 释放每次
			free(residuals); residuals = NULL;
			for (size_t i = 0; i < block_sizes.size(); i++)
			{
				free(jacobians[i]); jacobians[i] = NULL;
			}
			jacobians = NULL;
		}

		// 添加残差块
		ceres::ResidualBlockId blkId = problem_->AddResidualBlock(cost_function, NULL, vec);
		return blkId;
	}

	bool TrajectoryEstimator::AddLineMeasurements(const LiDMeasType type, const double dst_weight, const double ang_weight)
	{
		if (pFrameContainer == NULL || pFrameContainer->MapLines.size() <= 0)
			return false;

		// 遍历所有地图点
		int		nLine = 0, nObs = 0;
		std::map<int64_t, MapLine*>::iterator iter;
		for (iter = pFrameContainer->MapLines.begin(); iter != pFrameContainer->MapLines.end(); ++iter)
		{
			iter->second->residuals_blkid.clear();
			if (iter->second->m_pFrames.size() < 2)
				continue;

			std::vector<std::pair<lpostk::Frame*, size_t>> vec(iter->second->m_pFrames.begin(), iter->second->m_pFrames.end());
			std::sort(vec.begin(), vec.end(), CompareFrameIdPair);
			nLine++;

			if (type == LiDMeasType::LOCALTOGLOBAL) // 每一帧与全局系
			{
				for (size_t i = 0; i < vec.size(); i++)
				{
					size_t innerId = vec[i].second;
					LineObs obs = vec[i].first->Lines[innerId];
					obs.timestamp = fmod(obs.timestamp, 604800.0);
					ceres::ResidualBlockId blkId = AddLineLocalToGlobalMeasurement(obs, &iter->second->Alpha, &iter->second->Rotation, dst_weight, ang_weight);

					ObsPair	obspair;
					obspair.Id1 = innerId; obspair.pFrame1 = vec[i].first;
					iter->second->residuals_blkid.push_back(std::pair<ceres::ResidualBlockId, ObsPair>(blkId, obspair));
					nObs++;
				}
			}
			else if (type == LiDMeasType::LOCALTOLOCAL) // 任意两帧之间
			{
				for (size_t i = 0; i < vec.size() - 1; i++)
				{
					// 第一帧
					size_t innerId1 = vec[i].second;
					LineObs obs1 = vec[i].first->Lines[innerId1];

					// 第二帧
					for (size_t j = i + 1; j < vec.size(); j++)
					{
						size_t innerId2 = vec[j].second;
						LineObs obs2 = vec[j].first->Lines[innerId2];

						// 只保留周秒
						obs1.timestamp = fmod(obs1.timestamp, 604800.0);
						obs2.timestamp = fmod(obs2.timestamp, 604800.0);

						//if (j - i < 20)	continue;
						ceres::ResidualBlockId blkId = AddLineLocalToLocalMeasurement(obs1, obs2, dst_weight, ang_weight);

						ObsPair	obspair;
						obspair.Id1 = innerId1; obspair.Id2 = innerId2;
						obspair.pFrame1 = vec[i].first; obspair.pFrame2 = vec[j].first;
						iter->second->residuals_blkid.push_back(std::pair<ceres::ResidualBlockId, ObsPair>(blkId, obspair));
						nObs++;
					}
				}
			}
		}

		fprintf(stdout, "AddLineMeasurements nLine = %d nObs = %d\n", nLine, nObs);
		return nLine > 0;
	}

	ceres::ResidualBlockId TrajectoryEstimator::AddLineLocalToGlobalMeasurement(const LineObs &obs, double *Alpha, liso::SO3d *Rotation, double dst_weight, const double ang_weight)
	{
		// 计算样条节点位置
		const double time = obs.timestamp;
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {time, time} }, spline_meta);

		LineLocalToGlobalFactor* functor = new LineLocalToGlobalFactor(obs, spline_meta, dst_weight, ang_weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<LineLocalToGlobalFactor>(functor);
		cost_function->SetNumResiduals(2);

		// 参数块大小
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 姿态
			cost_function->AddParameterBlock(4);
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 位置
			cost_function->AddParameterBlock(3);
		cost_function->AddParameterBlock(4); // 外参R
		cost_function->AddParameterBlock(3); // 外参T
		cost_function->AddParameterBlock(1); // 直线尺度参数
		cost_function->AddParameterBlock(4); // 直线旋转参数

		// 参数块列表
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);       // 姿态
		AddControlPoints(spline_meta, vec, true); // 位置
		vec.emplace_back(calib_param_->so3_LtoI.data()); problem_->AddParameterBlock(calib_param_->so3_LtoI.data(), 4, local_parameterization); // 外参R
		vec.emplace_back(calib_param_->p_LinI.data()); problem_->AddParameterBlock(calib_param_->p_LinI.data(), 3); // 外参T	
		vec.emplace_back(Alpha); problem_->AddParameterBlock(Alpha, 1); // 直线尺度参数
		vec.emplace_back(Rotation->data()); problem_->AddParameterBlock(Rotation->data(), 4, local_parameterization); // 直线旋转参数

		// 添加残差块
		ceres::ResidualBlockId BlockId = problem_->AddResidualBlock(cost_function, NULL, vec);
		return BlockId;
	}
	ceres::ResidualBlockId TrajectoryEstimator::AddLineLocalToLocalMeasurement(const LineObs &obs1, const LineObs &obs2, double dst_weight, const double ang_weight)
	{
		// 计算样条节点位置
		const double time1 = obs1.timestamp;
		const double time2 = obs2.timestamp;
		double	t1 = time1 <= time2 ? time1 : time2;
		double	t2 = time1 >= time2 ? time1 : time2;
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {t1, t2} }, spline_meta);

		LineLocalToLocalFactor* functor = new LineLocalToLocalFactor(obs1, obs2, spline_meta, dst_weight, ang_weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<LineLocalToLocalFactor>(functor);
		cost_function->SetNumResiduals(2);

		// 参数块大小
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 姿态
			cost_function->AddParameterBlock(4);
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 位置
			cost_function->AddParameterBlock(3);
		cost_function->AddParameterBlock(4); // 外参R
		cost_function->AddParameterBlock(3); // 外参T

		// 参数块列表
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);       // 第一帧姿态
		AddControlPoints(spline_meta, vec, true); // 第一帧位置
		vec.emplace_back(calib_param_->so3_LtoI.data()); problem_->AddParameterBlock(calib_param_->so3_LtoI.data(), 4, local_parameterization); // 外参R
		vec.emplace_back(calib_param_->p_LinI.data()); problem_->AddParameterBlock(calib_param_->p_LinI.data(), 3); // 外参T	

		// 添加残差块
		ceres::ResidualBlockId BlockId = problem_->AddResidualBlock(cost_function, NULL, vec);
		return BlockId;
	}

	void TrajectoryEstimator::printPoleResiduals()
	{

	}

	bool TrajectoryEstimator::AddPoleMeasurements(const LiDMeasType type, const double dst_weight, const double ang_weight)
	{
		if (pFrameContainer == NULL || pFrameContainer->MapPlanes.size() <= 0)
			return false;

		// 遍历所有地图点
		int		nPole = 0, nObs = 0;
		//std::map<int64_t, MapPlane*>::iterator iter;
		//for (iter = pFrameContainer->MapPlanes.begin(); iter != pFrameContainer->MapPlanes.end(); ++iter)
		//{
		//	iter->second->residuals_blkid.clear();
		//	if (iter->second->m_pFrames.size() < 2)
		//		continue;

		//	std::vector<std::pair<lpostk::Frame*, size_t>> vec(iter->second->m_pFrames.begin(), iter->second->m_pFrames.end());
		//	std::sort(vec.begin(), vec.end(), CompareFrameIdPair);
		//	nPlane++;

		//	if (type == LiDMeasType::LOCALTOGLOBAL) // 每一帧与全局系
		//	{
		//		for (size_t i = 0; i < vec.size(); i++)
		//		{
		//			size_t innerId = vec[i].second;
		//			PlaneObs obs = vec[i].first->Planes[innerId];
		//			obs.timestamp = fmod(obs.timestamp, 604800.0);
		//			ceres::ResidualBlockId blkId = AddPlaneLocalToGlobalMeasurement(obs, iter->second->CPW, dst_weight, ang_weight);

		//			ObsPair	obspair;
		//			obspair.Id1 = innerId; obspair.pFrame1 = vec[i].first;
		//			iter->second->residuals_blkid.push_back(std::pair<ceres::ResidualBlockId, ObsPair>(blkId, obspair));
		//			nObs++;
		//		}
		//	}
		//	else if (type == LiDMeasType::LOCALTOLOCAL) // 任意两帧之间
		//	{
		//		for (size_t i = 0; i < vec.size() - 1; i++)
		//		{
		//			// 第一帧
		//			size_t innerId1 = vec[i].second;
		//			PlaneObs obs1 = vec[i].first->Planes[innerId1];

		//			// 第二帧
		//			for (size_t j = i + 1; j < vec.size(); j++)
		//			{
		//				size_t innerId2 = vec[j].second;
		//				PlaneObs obs2 = vec[j].first->Planes[innerId2];

		//				// 只保留周秒
		//				obs1.timestamp = fmod(obs1.timestamp, 604800.0);
		//				obs2.timestamp = fmod(obs2.timestamp, 604800.0);

		//				//if (j - i < 20)	continue;
		//				ceres::ResidualBlockId blkId = AddPlaneLocalToLocalMeasurement(obs1, obs2, dst_weight, ang_weight);

		//				ObsPair	obspair;
		//				obspair.Id1 = innerId1; obspair.Id2 = innerId2;
		//				obspair.pFrame1 = vec[i].first; obspair.pFrame2 = vec[j].first;
		//				iter->second->residuals_blkid.push_back(std::pair<ceres::ResidualBlockId, ObsPair>(blkId, obspair));
		//				nObs++;
		//			}
		//		}
		//	}
		//}

		fprintf(stdout, "AddPoleMeasurements nPlane = %d nObs = %d\n", nPole, nObs);
		return nPole > 0;
	}

	ceres::ResidualBlockId TrajectoryEstimator::AddPoleLocalToGlobalMeasurement(const PoleObs &obs, double *Alpha, liso::SO3d *Rotation, double dst_weight, const double ang_weight)
	{
		// 计算样条节点位置
		const double time = obs.timestamp;
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {time, time} }, spline_meta);

		PoleLocalToGlobalFactor* functor = new PoleLocalToGlobalFactor(obs, spline_meta, dst_weight, ang_weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<PoleLocalToGlobalFactor>(functor);
		cost_function->SetNumResiduals(2); // 点到直线的距离

		// 参数块大小
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 姿态
			cost_function->AddParameterBlock(4);
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 位置
			cost_function->AddParameterBlock(3);
		cost_function->AddParameterBlock(4); // 外参R
		cost_function->AddParameterBlock(3); // 外参T
		cost_function->AddParameterBlock(1); // 轴线尺度参数
		cost_function->AddParameterBlock(4); // 轴线旋转参数

		// 参数块列表
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);       // 姿态
		AddControlPoints(spline_meta, vec, true); // 位置
		vec.emplace_back(calib_param_->so3_LtoI.data()); problem_->AddParameterBlock(calib_param_->so3_LtoI.data(), 4, local_parameterization); // 外参R
		vec.emplace_back(calib_param_->p_LinI.data()); problem_->AddParameterBlock(calib_param_->p_LinI.data(), 3); // 外参T	
		vec.emplace_back(Alpha); problem_->AddParameterBlock(Alpha, 1); // 轴线尺度参数
		vec.emplace_back(Rotation->data()); problem_->AddParameterBlock(Rotation->data(), 4, local_parameterization); // 轴线旋转参数

		// 添加残差块
		ceres::ResidualBlockId BlockId = problem_->AddResidualBlock(cost_function, NULL, vec);
		return BlockId;
	}

	ceres::ResidualBlockId TrajectoryEstimator::AddPoleLocalToLocalMeasurement(const PoleObs &obs1, const PoleObs &obs2, double dst_weight, const double ang_weight)
	{
		// 计算样条节点位置
		const double time1 = obs1.timestamp;
		const double time2 = obs2.timestamp;
		double	t1 = time1 <= time2 ? time1 : time2;
		double	t2 = time1 >= time2 ? time1 : time2;
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {t1, t2} }, spline_meta);

		PoleLocalToLocalFactor* functor = new PoleLocalToLocalFactor(obs1, obs2, spline_meta, dst_weight, ang_weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<PoleLocalToLocalFactor>(functor);
		cost_function->SetNumResiduals(2); // 2 还是 4 

		// 参数块大小
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 姿态
			cost_function->AddParameterBlock(4);
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 位置
			cost_function->AddParameterBlock(3);
		cost_function->AddParameterBlock(4); // 外参R
		cost_function->AddParameterBlock(3); // 外参T

		// 参数块列表
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);       // 姿态
		AddControlPoints(spline_meta, vec, true); // 位置
		vec.emplace_back(calib_param_->so3_LtoI.data()); problem_->AddParameterBlock(calib_param_->so3_LtoI.data(), 4, local_parameterization); // 外参R
		vec.emplace_back(calib_param_->p_LinI.data()); problem_->AddParameterBlock(calib_param_->p_LinI.data(), 3); // 外参T	

		// 添加残差块
		ceres::ResidualBlockId BlockId = problem_->AddResidualBlock(cost_function, NULL, vec);
		return BlockId;
	}

	void TrajectoryEstimator::printPlaneResiduals()
	{

	}

	bool TrajectoryEstimator::AddPlaneMeasurements(const LiDMeasType type, const double dst_weight, const double ang_weight)
	{
		if (pFrameContainer == NULL || pFrameContainer->MapPlanes.size() <= 0)
			return false;

		// 遍历所有地图点
		int		nPlane = 0, nObs = 0;
		std::map<int64_t, MapPlane*>::iterator iter;
		for (iter = pFrameContainer->MapPlanes.begin(); iter != pFrameContainer->MapPlanes.end(); ++iter)
		{
			iter->second->residuals_blkid.clear();
			if (iter->second->m_pFrames.size() < 2)
				continue;

			std::vector<std::pair<lpostk::Frame*, size_t>> vec(iter->second->m_pFrames.begin(), iter->second->m_pFrames.end());
			std::sort(vec.begin(), vec.end(), CompareFrameIdPair);
			nPlane++;

			if (type == LiDMeasType::LOCALTOGLOBAL) // 每一帧与全局系
			{
				for (size_t i = 0; i < vec.size(); i++)
				{
					size_t innerId = vec[i].second;
					PlaneObs obs = vec[i].first->Planes[innerId];
					obs.timestamp = fmod(obs.timestamp, 604800.0);
					ceres::ResidualBlockId blkId = AddPlaneLocalToGlobalMeasurement(obs, iter->second->CPW, dst_weight, ang_weight);

					ObsPair	obspair;
					obspair.Id1 = innerId; obspair.pFrame1 = vec[i].first;
					iter->second->residuals_blkid.push_back(std::pair<ceres::ResidualBlockId, ObsPair>(blkId, obspair));
					nObs++;
				}
			}
			else if (type == LiDMeasType::LOCALTOLOCAL) // 任意两帧之间
			{
				for (size_t i = 0; i < vec.size() - 1; i++)
				{
					// 第一帧
					size_t innerId1 = vec[i].second;
					PlaneObs obs1 = vec[i].first->Planes[innerId1];

					// 第二帧
					for (size_t j = i + 1; j < vec.size(); j++)
					{
						size_t innerId2 = vec[j].second;
						PlaneObs obs2 = vec[j].first->Planes[innerId2];

						// 只保留周秒
						obs1.timestamp = fmod(obs1.timestamp, 604800.0);
						obs2.timestamp = fmod(obs2.timestamp, 604800.0);

						//if (j - i < 20)	continue;
						ceres::ResidualBlockId blkId = AddPlaneLocalToLocalMeasurement(obs1, obs2, dst_weight, ang_weight);

						ObsPair	obspair;
						obspair.Id1 = innerId1; obspair.Id2 = innerId2;
						obspair.pFrame1 = vec[i].first; obspair.pFrame2 = vec[j].first;
						iter->second->residuals_blkid.push_back(std::pair<ceres::ResidualBlockId, ObsPair>(blkId, obspair));
						nObs++;
					}
				}
			}
		}

		fprintf(stdout, "AddPlaneMeasurements nPlane = %d nObs = %d\n", nPlane, nObs);
		return nPlane > 0;
	}

	ceres::ResidualBlockId TrajectoryEstimator::AddPlaneLocalToGlobalMeasurement(const PlaneObs &obs, double CPW[3], double dst_weight, const double ang_weight)
	{
		// 计算样条节点位置
		const double time = obs.timestamp;
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {time, time} }, spline_meta);

		PlaneLocalToGlobalCPFactor* functor = new PlaneLocalToGlobalCPFactor(obs, spline_meta, dst_weight, ang_weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<PlaneLocalToGlobalCPFactor>(functor);
		cost_function->SetNumResiduals(2); // 距离 + 角度

		// 参数块大小
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 当前帧姿态
			cost_function->AddParameterBlock(4);
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 当前帧位置
			cost_function->AddParameterBlock(3);
		cost_function->AddParameterBlock(4); // 外参R
		cost_function->AddParameterBlock(3); // 外参T
		cost_function->AddParameterBlock(3); // 平面最近点

		// 参数块列表
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);       // 当前帧姿态
		AddControlPoints(spline_meta, vec, true); // 当前帧位置
		vec.emplace_back(calib_param_->so3_LtoI.data()); problem_->AddParameterBlock(calib_param_->so3_LtoI.data(), 4, local_parameterization); // 外参R
		vec.emplace_back(calib_param_->p_LinI.data()); problem_->AddParameterBlock(calib_param_->p_LinI.data(), 3); // 外参T	
		vec.emplace_back(CPW); // 平面最近点

		// 添加残差块
		ceres::ResidualBlockId BlockId = problem_->AddResidualBlock(cost_function, NULL, vec);
		return BlockId;
	}

	ceres::ResidualBlockId TrajectoryEstimator::AddPlaneLocalToLocalMeasurement(const PlaneObs &obs1, const PlaneObs &obs2, double dst_weight, const double ang_weight)
	{
		// 计算样条节点位置
		const double time1 = obs1.timestamp;
		const double time2 = obs2.timestamp;
		double	t1 = time1 <= time2 ? time1 : time2;
		double	t2 = time1 >= time2 ? time1 : time2;
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {t1, t2} }, spline_meta);

		PlaneLocalToLocalFactor* functor = new PlaneLocalToLocalFactor(obs1, obs2, spline_meta, dst_weight, ang_weight);
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<PlaneLocalToLocalFactor>(functor);
		cost_function->SetNumResiduals(4); // 

		// 参数块大小
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 姿态
			cost_function->AddParameterBlock(4);
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 位置
			cost_function->AddParameterBlock(3);
		cost_function->AddParameterBlock(4); // 外参R
		cost_function->AddParameterBlock(3); // 外参T

		// 参数块列表
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);       // 姿态
		AddControlPoints(spline_meta, vec, true); // 位置
		vec.emplace_back(calib_param_->so3_LtoI.data()); problem_->AddParameterBlock(calib_param_->so3_LtoI.data(), 4, local_parameterization); // 外参R
		vec.emplace_back(calib_param_->p_LinI.data()); problem_->AddParameterBlock(calib_param_->p_LinI.data(), 3); // 外参T	

		// 添加残差块
		ceres::ResidualBlockId BlockId = problem_->AddResidualBlock(cost_function, NULL, vec);
		return BlockId;
	}


	/// @brief GNSS位置/速度约束
	///
	/// @param[in] obs			GNSS解算结果
	/// @param[in] blarm		GNSS杆臂
	/// @param[in] weight		先验位置/速度的权重
	/// @return void
	ceres::ResidualBlockId TrajectoryEstimator::addGNSSPosMeasurement(const GNSSPVA &obs, const Eigen::Vector3d &blarm, double weight)
	{
		// 计算样条节点位置
		const double time1 = obs.timestamp;
		const double time2 = obs.timestamp;
		double	t1 = time1 <= time2 ? time1 : time2;
		double	t2 = time1 >= time2 ? time1 : time2;
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {t1, t2} }, spline_meta);

		GNSSPositionFactor* functor = new GNSSPositionFactor(obs.timestamp, obs.position, blarm, spline_meta, Eigen::Vector3d(obs.poscov(0,0), obs.poscov(1, 1), obs.poscov(2, 2)));
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<GNSSPositionFactor>(functor);
		cost_function->SetNumResiduals(3); // 

		// 参数块大小
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 姿态
			cost_function->AddParameterBlock(4);
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 位置
			cost_function->AddParameterBlock(3);

		// 参数块列表
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);       // 姿态
		AddControlPoints(spline_meta, vec, true); // 位置

		// 添加残差块
		ceres::ResidualBlockId BlockId = problem_->AddResidualBlock(cost_function, NULL, vec);
		return BlockId;
	}


	ceres::ResidualBlockId TrajectoryEstimator::addGNSSVelMeasurement(const GNSSPVA &obs, const Eigen::Vector3d &blarm, double weight)
	{
		// 计算样条节点位置
		const double time1 = obs.timestamp;
		const double time2 = obs.timestamp;
		double	t1 = time1 <= time2 ? time1 : time2;
		double	t2 = time1 >= time2 ? time1 : time2;
		SplineMeta<SplineOrder> spline_meta;
		trajectory_->CaculateSplineMeta({ {t1, t2} }, spline_meta);

		GNSSVelocityFactor* functor = new GNSSVelocityFactor(obs.timestamp, obs.velocity, blarm, spline_meta, Eigen::Vector3d(obs.velcov(0, 0), obs.velcov(1, 1), obs.velcov(2, 2)));
		auto* cost_function = new ceres::DynamicAutoDiffCostFunction<GNSSVelocityFactor>(functor);
		cost_function->SetNumResiduals(3); // 

		// 参数块大小
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 姿态
			cost_function->AddParameterBlock(4);
		for (int i = 0; i < spline_meta.NumParameters(); i++) // 位置
			cost_function->AddParameterBlock(3);

		// 参数块列表
		std::vector<double*> vec;
		AddControlPoints(spline_meta, vec);       // 姿态
		AddControlPoints(spline_meta, vec, true); // 位置

		// 添加残差块
		ceres::ResidualBlockId BlockId = problem_->AddResidualBlock(cost_function, NULL, vec);
		return BlockId;

	}


	void TrajectoryEstimator::SetTrajectorControlPointVariable(double min_time,
		double max_time) {
		size_t min_s = trajectory_->computeTIndex(min_time).second;
		size_t max_s = trajectory_->computeTIndex(max_time).second;

		std::cout << "[SetTrajectorControlPointVariable]: " << min_s << ", "
			<< max_s + SplineOrder - 1 << "\n";
		for (size_t i = min_s; i < (max_s + SplineOrder); i++) {
			problem_->AddParameterBlock(trajectory_->getKnotSO3(i).data(), 4,
				local_parameterization);
			problem_->SetParameterBlockVariable(trajectory_->getKnotSO3(i).data());

			problem_->AddParameterBlock(trajectory_->getKnotPos(i).data(), 3);
			problem_->SetParameterBlockVariable(trajectory_->getKnotPos(i).data());
		}
	}

	/// Add callback for debug
	void TrajectoryEstimator::AddCallback(bool needs_state, std::string filename) {
		std::unique_ptr<CheckStateCallback> cb =
			std::make_unique<CheckStateCallback>(filename);

		cb->addCheckState("px py pz", 3, calib_param_->p_LinI.data());
		cb->addCheckState("qx qy qz qw", 4, calib_param_->so3_LtoI.data());

		for (SegmentCalibParam& v : calib_param_->segment_param) {
			if (problem_->HasParameterBlock(v.acce_bias.data())) {
				cb->addCheckState("ax ay az", 3, v.acce_bias.data());
			}
			if (problem_->HasParameterBlock(v.gyro_bias.data())) {
				cb->addCheckState("gx gy gz", 3, v.gyro_bias.data());
			}
			if (problem_->HasParameterBlock(v.g_refine.data())) {
				cb->addCheckState("g_refinex g_refiney", 2, v.g_refine.data());
			}
			if (problem_->HasParameterBlock(&v.time_offset)) {
				cb->addCheckState("time_offset", 1, &v.time_offset);
			}
		}
		callbacks_.push_back(std::move(cb));

		// If any callback requires state, the flag must be set
		if (!callback_needs_state_) callback_needs_state_ = needs_state;
	}

	/// @brief 锁定轨迹节点
	///
	/// @param[in] lock_P    锁定位置
	/// @param[in] lock_R    锁定姿态
	/// @return void
	void TrajectoryEstimator::LockTrajKnots(bool lock_P, bool lock_R)
	{
		if (!lock_P && !lock_R)
			return;

		for (size_t i = 0; i < trajectory_->numKnots(); i++)
		{
			// 位姿
			if (lock_P && problem_->HasParameterBlock(trajectory_->getKnotPos(i).data()))
				problem_->SetParameterBlockConstant(trajectory_->getKnotPos(i).data());

			// 姿态
			if (lock_R && problem_->HasParameterBlock(trajectory_->getKnotSO3(i).data()))
				problem_->SetParameterBlockConstant(trajectory_->getKnotSO3(i).data());
		}

		return;
	}

	/// @brief 锁定IMU参数
	///
	/// @param[in] lock_ab    锁定加计零偏
	/// @param[in] lock_wb    锁定陀螺零偏
	/// @param[in] lock_g     锁定重力
	/// @return void
	void TrajectoryEstimator::LockIMUState(bool lock_ab, bool lock_wb, bool lock_g)
	{
		if (calib_param_ == nullptr || calib_param_->segment_param.size() == 0)
			return;

		for (SegmentCalibParam& v : calib_param_->segment_param)
		{
			// 加计零偏
			if (lock_ab && problem_->HasParameterBlock(v.acce_bias.data()))
			{
				//problem_->AddParameterBlock(v.acce_bias.data(), 3);
				problem_->SetParameterBlockConstant(v.acce_bias.data());
			}

			// 陀螺零偏
			if (lock_wb && problem_->HasParameterBlock(v.gyro_bias.data()))
			{
				//problem_->AddParameterBlock(v.gyro_bias.data(), 3);
				problem_->SetParameterBlockConstant(v.gyro_bias.data());
			}

			// 重力方向
			if (lock_g && problem_->HasParameterBlock(v.g_refine.data())) {
				//problem_->AddParameterBlock(v.g_refine.data(), 2);
				problem_->SetParameterBlockConstant(v.g_refine.data());
			}
		}
	}

	/// @brief 锁定激光惯导外参
	///
	/// @param[in] lock_P    锁定杆臂
	/// @param[in] lock_R    锁定旋转
	/// @return void
	void TrajectoryEstimator::LockExtrinsicParam(bool lock_P, bool lock_R)
	{
		// 杆臂
		if (lock_P && problem_->HasParameterBlock(calib_param_->p_LinI.data()))
		{
			//problem_->AddParameterBlock(calib_param_->p_LinI.data(), 3);
			problem_->SetParameterBlockConstant(calib_param_->p_LinI.data());
		}

		// 旋转
		if (lock_R && problem_->HasParameterBlock(calib_param_->so3_LtoI.data()))
		{
			//problem_->AddParameterBlock(calib_param_->so3_LtoI.data(), 4, local_parameterization);
			problem_->SetParameterBlockConstant(calib_param_->so3_LtoI.data());
		}
	}

	void TrajectoryEstimator::ShowParameters()
	{
		const bool showKnots = true;
		bool	bHasPos = false, bFixPos = false, bHasSO3 = false, bFixSO3 = false;
		bool	bHasAcc = false, bFixAcc = false, bHasGyr = false, bFixGyr = false, bHasGra = false, bFixGra = false;

		// 轨迹节点
		if (showKnots)
		{
			for (size_t i = 0; i < trajectory_->numKnots(); i++)
			{
				bHasPos = problem_->HasParameterBlock(trajectory_->getKnotPos(i).data());
				bFixPos = bHasPos ? problem_->IsParameterBlockConstant(trajectory_->getKnotPos(i).data()) : false;
				bHasSO3 = problem_->HasParameterBlock(trajectory_->getKnotSO3(i).data());
				bFixSO3 = bHasSO3 ? problem_->IsParameterBlockConstant(trajectory_->getKnotSO3(i).data()) : false;
				fprintf(stdout, "Knot %10zd %d %d %d %d\n", i, bHasPos, bFixPos, bHasSO3, bFixSO3);
			}
		}

		// 分段惯导参数
		if (calib_param_)
		{
			for (size_t i = 0; i < calib_param_->segment_param.size(); i++)
			{
				SegmentCalibParam& v = calib_param_->segment_param[i];

				// 陀螺零偏
				bHasGyr = problem_->HasParameterBlock(v.gyro_bias.data());
				bFixGyr = bHasGyr ? problem_->IsParameterBlockConstant(v.gyro_bias.data()) : false;

				// 加计零偏
				bHasAcc = problem_->HasParameterBlock(v.acce_bias.data());
				bFixAcc = bHasAcc ? problem_->IsParameterBlockConstant(v.acce_bias.data()) : false;

				// 重力方向
				bHasGra = problem_->HasParameterBlock(v.g_refine.data());
				bFixGra = bHasGra ? problem_->IsParameterBlockConstant(v.g_refine.data()) : false;

				fprintf(stdout, "Segment %10zd %d %d %d %d %d %d\n", i, bHasGyr, bFixGyr, bHasAcc, bFixAcc, bHasGra, bFixGra);
			}
		}

		// 激光惯导外参
		bHasPos = problem_->HasParameterBlock(calib_param_->p_LinI.data());
		bFixPos = bHasPos ? problem_->IsParameterBlockConstant(calib_param_->p_LinI.data()) : false;
		bHasSO3 = problem_->HasParameterBlock(calib_param_->so3_LtoI.data());
		bFixSO3 = bHasPos ? problem_->IsParameterBlockConstant(calib_param_->so3_LtoI.data()) : false;
		fprintf(stdout, "LiDAR Extrinsic %d %d %d %d\n", bHasPos, bFixPos, bHasSO3, bFixSO3);

		// 地图点
		bool	bHas = false, bFix = false;
		if (pFrameContainer)
		{
			// 点
			std::map<int64_t, MapPoint*>::iterator iter_point;
			for (iter_point = pFrameContainer->MapPoints.begin(); iter_point != pFrameContainer->MapPoints.end(); iter_point++)
			{
				bHas = problem_->HasParameterBlock(iter_point->second->XYZ);
				bFix = bHas ? problem_->IsParameterBlockConstant(iter_point->second->XYZ) : false;
				fprintf(stdout, "MapPoint %010lld %d %d\n", iter_point->first, bHas, bFix);
			}

			// 线
			std::map<int64_t, MapLine*>::iterator iter_line;
			for (iter_line = pFrameContainer->MapLines.begin(); iter_line != pFrameContainer->MapLines.end(); iter_line++)
			{
				bHas = problem_->HasParameterBlock(&iter_line->second->Alpha);
				bFix = bHas ? problem_->IsParameterBlockConstant(&iter_line->second->Alpha) : false;
				bHasSO3 = problem_->HasParameterBlock(iter_line->second->Rotation.data());
				bFixSO3 = bHasSO3 ? problem_->IsParameterBlockConstant(iter_line->second->Rotation.data()) : false;
				fprintf(stdout, "MapLine  %010lld %d %d %d %d\n", iter_point->first, bHas, bFix, bHasSO3, bFixSO3);
			}

			// 柱
			std::map<int64_t, MapPole*>::iterator iter_pole;
			for (iter_pole = pFrameContainer->MapPoles.begin(); iter_pole != pFrameContainer->MapPoles.end(); iter_pole++)
			{
				bHas = problem_->HasParameterBlock(&iter_pole->second->Alpha);
				bFix = bHas ? problem_->IsParameterBlockConstant(&iter_pole->second->Alpha) : false;
				bHasSO3 = problem_->HasParameterBlock(iter_pole->second->Rotation.data());
				bFixSO3 = bHasSO3 ? problem_->IsParameterBlockConstant(iter_pole->second->Rotation.data()) : false;
				fprintf(stdout, "MapPole  %010lld %d %d %d %d\n", iter_point->first, bHas, bFix, bHasSO3, bFixSO3);
			}

			// 面
			std::map<int64_t, MapPlane*>::iterator iter_plane;
			for (iter_plane = pFrameContainer->MapPlanes.begin(); iter_plane != pFrameContainer->MapPlanes.end(); iter_plane++)
			{
				bHas = problem_->HasParameterBlock(iter_plane->second->CPW);
				bFix = bHas ? problem_->IsParameterBlockConstant(iter_plane->second->CPW) : false;
				fprintf(stdout, "MapPlane %010lld %d %d\n", iter_point->first, bHas, bFix);
			}
		}

		return;
	}

	ceres::Solver::Summary TrajectoryEstimator::Solve(int max_iterations, bool progress, int num_threads)
	{
		// 锁定轨迹节点
		LockTrajKnots(options_.lock_trajP, options_.lock_trajR);

		// 锁定零偏/陀螺零偏/重力
		LockIMUState(options_.lock_ab, options_.lock_wb, options_.local_traj || options_.lock_g);

		// 锁定激光惯导外参
		LockExtrinsicParam(options_.lock_tlb, options_.lock_Rlb);

		// 输出参数的估计状态
		ShowParameters();

		printf("ParameterBlocks = %d\n", problem_->NumParameterBlocks());
		printf("ResidualBlocks  = %d\n", problem_->NumResidualBlocks());

		ceres::Solver::Options options;

		options.minimizer_type = ceres::TRUST_REGION;
		//  options.gradient_tolerance = 0.01 * Sophus::Constants<double>::epsilon();
		//  options.function_tolerance = 0.01 * Sophus::Constants<double>::epsilon();
		options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
		options.linear_solver_type = ceres::DENSE_NORMAL_CHOLESKY;
		//    options.trust_region_strategy_type = ceres::DOGLEG;
		//    options.dogleg_type = ceres::SUBSPACE_DOGLEG;

		//    options.linear_solver_type = ceres::SPARSE_SCHUR;
		options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
		options.minimizer_progress_to_stdout = true;

		if (!options_.lock_tlb) {
			// function_tolerance = 1e-6;
			// gradient_tolerance = 1e-10;
			// parameter_tolerance = 1e-8;
			// options.function_tolerance = 1e-20;
			// options.parameter_tolerance = 1e-20;
			// options.initial_trust_region_radius = 1e3;
			max_iterations = 100;
		}

		options.minimizer_progress_to_stdout = progress;

		if (num_threads < 1) {
			num_threads = std::thread::hardware_concurrency();
		}
		options.num_threads = num_threads;
		options.max_num_iterations = max_iterations;

		if (callbacks_.size() > 0) {
			for (auto& cb : callbacks_) {
				options.callbacks.push_back(cb.get());
			}

			if (callback_needs_state_) options.update_state_every_iteration = true;
		}

		if (1)
		{
			options.minimizer_type = ceres::TRUST_REGION;
			options.linear_solver_type = ceres::SPARSE_SCHUR;
			//options.max_num_iterations = 1000;
			options.gradient_tolerance = 1e-14;
			options.function_tolerance = 1e-14;
			options.parameter_tolerance = 1.0e-14;
			options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
			options.update_state_every_iteration = true;
			options.minimizer_progress_to_stdout = true;
		}

		ceres::Solver::Summary summary;
		ceres::Solve(options, problem_.get(), &summary);
		//std::cout << summary.BriefReport() << std::endl;

		// update state
		if (calib_param_)
		{
			calib_param_->UpdateExtrinicParam();
			calib_param_->UpdateGravity();
		}
		//  getCovariance();
		return summary;
	}

	bool TrajectoryEstimator::getCovariance() {
		ceres::Covariance::Options options;
		// options.algorithm_type = ceres::DENSE_SVD;
		options.apply_loss_function = false;
		ceres::Covariance covariance(options);

		if (!problem_->HasParameterBlock(calib_param_->p_LinI.data()) ||
			!problem_->HasParameterBlock(calib_param_->so3_LtoI.data())) {
			return false;
		}

		std::vector<const double*> vec;
		vec.emplace_back(calib_param_->p_LinI.data());
		vec.emplace_back(calib_param_->so3_LtoI.data());

		if (!covariance.Compute(vec, problem_.get())) {
			std::cout << "[CovarianceMatrixInTangentSpace] J^TJ is a rank deficient matrix\n";
			return false;
		}

		double m2cm = 100;
		double rad2degree = 180.0 / M_PI;

		Eigen::Matrix<double, 6, 6> cov = Eigen::Matrix<double, 6, 6>::Zero();
		covariance.GetCovarianceMatrixInTangentSpace(vec, cov.data());

		Eigen::VectorXd diag;
		diag = cov.diagonal();
		diag.head<6>() = diag.head<6>().cwiseSqrt();
		diag.segment<3>(0) *= m2cm;
		diag.segment<3>(3) *= rad2degree;

		std::cout << std::fixed << std::setprecision(9);
		std::cout << "[CovarianceMatrixInTangentSpace] \n" << cov << std::endl;
		std::cout << YELLOW;
		std::cout << "[std] pos (cm)    : " << diag.segment<3>(0).transpose()
			<< std::endl;
		std::cout << "[std] rot (degree): " << diag.segment<3>(3).transpose() << RESET
			<< std::endl;
		return true;
	}

}  // namespace liso
