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

#pragma once

#include <Eigen/Core>
#include <sophus/se3.hpp>
#include <sophus/so3.hpp>

namespace liso {

	const double GRAVITY_NORM = -9.797;

	namespace gravity_factor {

		/// @brief 更新重力
		///
		/// @param[in] g_param		俯仰角和横滚角
		/// @return 对应的惯导参数
		template <typename T>
		Eigen::Matrix<T, 3, 1> refined_gravity(Eigen::Map<const Eigen::Matrix<T, 2, 1>>& g_param)
		{
			T cr = ceres::cos(g_param[0]); // 横滚角
			T sr = ceres::sin(g_param[0]);
			T cp = ceres::cos(g_param[1]); // 俯仰角
			T sp = ceres::sin(g_param[1]);
			return Eigen::Matrix<T, 3, 1>(-sp * cr * T(GRAVITY_NORM),
				sr * T(GRAVITY_NORM),
				-cr * cp * T(GRAVITY_NORM));
		}
	}  // namespace gravity_factor

	using SO3d = Sophus::SO3<double>;
	using SE3d = Sophus::SE3<double>;

	struct IMUData {
		double timestamp;
		Eigen::Matrix<double, 3, 1> gyro;
		Eigen::Matrix<double, 3, 1> accel;
		Eigen::Quaterniond orientation;
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};

	struct PoseData {
		PoseData()
			: week(-1), sow(0), timestamp(0),
			position(Eigen::Vector3d(0, 0, 0)),
			orientation(SO3d(Eigen::Quaterniond::Identity())) {}

		int		week;	// 周
		double	sow;	// 周秒
		double timestamp;
		Eigen::Vector3d position;
		Eigen::Vector3d velocity;
		SO3d orientation; // Rbw/Rbe
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};

	struct OdomData {
		double timestamp;
		Eigen::Matrix4d pose;
	};

	inline bool	SlerpPoseData(const PoseData pose1, const PoseData pose2, const double t, PoseData &pose)
	{
		double	t1 = pose1.timestamp;
		double	t2 = pose2.timestamp;
		double	ratio = (t - t1) / (t2 - t1); // 这里可以看出当t1和t2过于接近时会引起数值误差

		Eigen::Vector3d		tran1 = pose1.position;
		Eigen::Quaterniond	quat1 = pose1.orientation.unit_quaternion();
		Eigen::Vector3d		tran2 = pose2.position;
		Eigen::Quaterniond	quat2 = pose2.orientation.unit_quaternion();

		Eigen::Vector3d		tran;
		Eigen::Quaterniond	quat;
		assert(t1 <= t && t <= t2);
		// 当t1和t2过于接近时，选择与当前时刻最近的位姿作为插值结果即可
		if (t2 - t1 <= std::numeric_limits<double>::epsilon()) {
			quat = ((t1 + t2) > 2 * t) ? quat1 : quat2;
			tran = ((t1 + t2) > 2 * t) ? tran1 : tran2;
		}
		else
		{
			quat = quat1.slerp(ratio, quat2);		// 旋转部分位姿使用球面线性插值
			tran = tran1 + ratio * (tran2 - tran1); // 平移部分位姿使用线性插值
		}

		pose.position = tran;
		pose.orientation = pose.orientation = liso::SO3d(quat);
		return true;
	}

}  // namespace liso
