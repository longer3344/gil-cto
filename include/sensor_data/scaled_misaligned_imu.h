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

#ifndef SCALED_MISALIGNED_IMU_H
#define SCALED_MISALIGNED_IMU_H

#include <ceres/ceres.h>
#include <Eigen/Eigen>
#include <Eigen/Geometry>

namespace liso {
	struct IMUIntrinsic {

		/// @brief 构造函数
		///
		///
		IMUIntrinsic()
		{
			// 陀螺尺度因子(0-2)和非正交误差(3-5)
			Ma_vec_.head<3>() = Eigen::Vector3d::Ones(); // 比例因子, 默认1.0
			Ma_vec_.tail<3>() = Eigen::Vector3d::Zero(); // 非正交误差, 默认0.0

			// 加计尺度因子(0-2)和非正交误差(3-5)
			Mw_vec_.head<3>() = Eigen::Vector3d::Ones(); // 比例因子, 默认1.0
			Mw_vec_.tail<3>() = Eigen::Vector3d::Zero(); // 非正交误差, 默认0.0

			Aw_vec_ = Eigen::Matrix<double, 9, 1>::Zero(); // 加计对陀螺的影响, 默认0.0
			q_WtoA_ = Eigen::Quaterniond::Identity();      // 陀螺相对于加计的安置角, 默认单位阵
		}

		/// @brief IMU原始数据修正
		///
		/// @param[in] Mw           陀螺比例因子和非正交误差
		/// @param[in] Aw           加计对陀螺的影响,一般为零矩阵
		/// @param[in] Ma           加计比例因子和非正交误差
		/// @param[in] q_WtoA       陀螺相对于加计的安置角,一般为单位阵
		/// @return void
		template <typename T>
		static void GetCalibrated(const Eigen::Matrix<T, 3, 3>& Mw, const Eigen::Matrix<T, 3, 3>& Aw, 
			const Eigen::Matrix<T, 3, 3>& Ma, const Eigen::Quaternion<T>& q_WtoA,
			Eigen::Matrix<T, 3, 1>& w_b, Eigen::Matrix<T, 3, 1>& a_b) 
		{
			w_b = (Mw * (q_WtoA * w_b) + Aw * a_b).eval();
			a_b = (Ma * a_b).eval();
		}

		/// @brief 获取陀螺比例因子和非正交误差
		///
		/// @return double[9]
		double* GetMwVector() { return Mw_vec_.data(); }

		/// @brief 获取加计比例因子和非正交误差
		///
		/// @return double[9]
		double* GetMaVector() { return Ma_vec_.data(); }

		/// @brief 获取加计对陀螺的影响
		///
		/// @return double[9]
		double* GetAwVector() { return Aw_vec_.data(); }

		/// @brief 获取陀螺相对于加计的安置角
		///
		/// @return double[9]
		double* GetQWtoAVector() { return q_WtoA_.coeffs().data(); }

		void ShowIMUParam() const 
		{
			std::cout << std::fixed << std::setprecision(6);
			std::cout << "Sw," << Mw_vec_[0] << "," << Mw_vec_[1] << "," << Mw_vec_[2] << std::endl; // 陀螺尺度因子
			std::cout << "Mw," << Mw_vec_[3] << "," << Mw_vec_[4] << "," << Mw_vec_[5] << std::endl; // 陀螺非正交误差
			std::cout << "Sa," << Ma_vec_[0] << "," << Ma_vec_[1] << "," << Ma_vec_[2] << std::endl; // 加计尺度因子
			std::cout << "Ma," << Ma_vec_[3] << "," << Ma_vec_[4] << "," << Ma_vec_[5] << std::endl; // 加计非正交误差

			// 加计对陀螺的影响
			std::cout << "Aw," << Aw_vec_[0];
			for (int i = 1; i < 9; ++i) std::cout << "," << Aw_vec_[i];
			std::cout << std::endl;

			// 陀螺相对于加计的安置角
			std::cout.unsetf(std::ios::fixed);
			Eigen::Vector3d euler_WtoA = (q_WtoA_.toRotationMatrix().eulerAngles(0, 1, 2)) * 180 / M_PI;
			std::cout << "euler_WtoA: " << euler_WtoA.transpose() << std::endl;
			Eigen::Vector3d euler_AtoW = (q_WtoA_.inverse().toRotationMatrix().eulerAngles(0, 1, 2)) * 180 / M_PI;
			std::cout << "euler_AtoW: " << euler_AtoW.transpose() << std::endl;
		}

	public:
		Eigen::Matrix<double, 6, 1> Mw_vec_; ///< 陀螺尺度因子(0-2)和非正交误差(3-5), 非正交矩阵是对称的,只需要3个值
		Eigen::Matrix<double, 6, 1> Ma_vec_; ///< 加计尺度因子(0-2)和非正交误差(3-5), 非正交矩阵是对称的,只需要3个值
		Eigen::Matrix<double, 9, 1> Aw_vec_; ///< 加计对陀螺的影响,一般为零矩阵
		Eigen::Quaterniond q_WtoA_;          ///< 陀螺相对于加计的安置角,一般为单位阵
	};

}  // namespace liso

#endif
