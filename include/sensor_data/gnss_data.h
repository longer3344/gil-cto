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


#ifndef GNSS_DATA_H
#define GNSS_DATA_H

#include <Eigen/Core>
#include <sophus/se3.hpp>
#include <sophus/so3.hpp>

namespace liso {

	struct GNSSPVA {
		GNSSPVA()
			: week(-1), sow(0), timestamp(0),
			position(Eigen::Vector3d(0, 0, 0)) {}

		int		week;	// 周
		double	sow;	// 周秒
		double	timestamp;
		Eigen::Vector3d position;
		Eigen::Matrix3d poscov;
		Eigen::Vector3d velocity;
		Eigen::Matrix3d velcov;
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};

	typedef struct tagELCPOSVEL
	{
		//GPSTIME gt;                        ///< GPS时间
		double  Pos[3], Vel[3], Att[3];    ///< 位置速度姿态 
		double  PosP[9], VelP[9], AttP[9]; ///< Covariance
		double  AttBL[3];                  ///< 双天线测姿的基线向量
		char    Q;                         ///< GNSS的Q因子
		char    AmbFix;                    ///< 模糊度固定情况
		float   DDOP;                      ///< 差分DDOP值

	} ELC_PVA;

}  // namespace liso

#endif
