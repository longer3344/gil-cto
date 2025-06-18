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

namespace liso {

struct TrajectoryEstimatorOptions {
  // If we should optimize the trajectory
  //bool lock_traj = false;

  // lock the extrinsic position/rotation/t_offset between imu and sensor
  bool lock_tlb = true;      ///< 是否不优化激光雷达杆臂
  bool lock_Rlb = true;      ///< 是否不优化激光雷达安置角
  bool lock_t_offset = true; ///< 是否不优化激光雷达时间偏差

  // If estimating the time offset, the max/min value of time offset
  double t_offset_padding = 0.02; ///< 激光雷达时间偏差阈值

  // lock the imu bias/gravity
  bool lock_ab = true; ///< 是否不优化加计零偏
  bool lock_wb = true; ///< 是否不优化陀螺零偏
  bool lock_g = true;  ///< 是否不优化重力
  bool lock_acc_scale = true;    ///< 是否不优化加计比例因子
  bool lock_acc_misalign = true; ///< 是否不优化加计非正交误差
  bool lock_gyr_scale = true;    ///< 是否不优化陀螺比例因子
  bool lock_gyr_misalign = true; ///< 是否不优化陀螺非正交误差
  bool lock_Aw = true;           ///< 是否不优化加计对陀螺的影响
  bool lock_R_WtoA = true;       ///< 是否不优化陀螺相对于加计的安置角

  bool lock_LiDAR_intrinsic = true; ///< 是否不优化雷达内参
  bool lock_IMU_intrinsic = true;   ///< 是否不优化IMU内参

  bool local_traj = true;   ///< 是否为局部轨迹,局部轨迹需要优化重力
  bool lock_trajP = false;	///< 是否不优化位置节点
  bool lock_trajR = false;  ///< 是否不优化姿态节点
  Eigen::Vector3d gravity = Eigen::Vector3d(0.0, 0.0, 0.0); ///< 重力加速度
};

}  // namespace liso
