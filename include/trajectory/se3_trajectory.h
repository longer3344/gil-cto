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

#ifndef SE3_TRAJECTORY_H
#define SE3_TRAJECTORY_H

#include <basalt/spline/se3_spline.h>
#include <sensor_data/calibration.h>
#include <sensor_data/cloud_type.h>
#include <sensor_data/gnss_data.h>
#include <sophus/se3.hpp>
#include <string>

namespace liso {

	class Trajectory : public basalt::Se3Spline<SplineOrder, double> {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
		typedef std::shared_ptr<Trajectory> Ptr;  ///< 智能指针

		// time_interval - 节点的时间间隔 start_time - 开始时间 segment_id - 轨迹分段ID
		Trajectory(double time_interval, double start_time = 0, size_t segment_id = 0)
			: basalt::Se3Spline<SplineOrder, double>(time_interval, start_time), segment_id_(segment_id)
		{
			time_interval_ = time_interval;

			// 开头插入N个空的数据体,最终节点的个数为N+P,N为阶数(次数+1),P为控制点个数
			this->extendKnotsTo(start_time, SO3d(Eigen::Quaterniond::Identity()), Eigen::Vector3d(0, 0, 0));
		}

		void SetCalibParam(std::shared_ptr<CalibParamManager> &calib_param) {
			//    calib_param_ = std::move(calib_param);
			calib_param_ = calib_param;
		}

		const std::shared_ptr<CalibParamManager> GetCalibParam() const {
			return calib_param_;
		}
		
		/// @brief 获取轨迹分段编号
		///
		/// @return 轨迹分段编号
		const size_t SegmentID() const { return segment_id_; }

		/// @brief 获取轨迹分段编号
		///
		/// @param[in] timestamp    时间
		/// @return 轨迹分段编号
		const size_t SegmentID(const double timestamp) const
		{
			// 只有一段轨迹
			if (calib_param_->segment_param.size() == 1 || calib_param_->segment_interval == 0.0
				|| timestamp <= calib_param_->traj_starttime)
			{
				return 0;
			}

			// 时间超过了
			if (timestamp >= calib_param_->traj_endtime)
			{
				return calib_param_->segment_param.size() - 1;
			}

			// 根据时间计算惯导参数索引
			double	dt = timestamp - calib_param_->traj_starttime;
			int		segment_id = (int)(dt / calib_param_->segment_interval);
			return segment_id;

			return segment_id_;
		}

		/// @brief 获取惯导相关参数
		///
		/// @return 对应的惯导参数
		SegmentCalibParam *GetTrajParam() {
			if (calib_param_->segment_param.size() == 1 || calib_param_->segment_interval == 0.0)
				return &calib_param_->segment_param[0];

			return &calib_param_->segment_param[segment_id_];
		}

		/// @brief 获取惯导相关参数
		///
		/// @param[in] timestamp    时间
		/// @return 对应的惯导参数
		SegmentCalibParam *GetTrajParam(const double timestamp) {
			int	segment_id = SegmentID(timestamp);
			return &calib_param_->segment_param[segment_id];
		}

		/// @brief 获取惯导相关参数
		///
		/// @return 对应的惯导参数
		const SegmentCalibParam *GetTrajParam() const {
			return &calib_param_->segment_param[segment_id_];
		}

		/// @brief 获取惯导相关参数
		///
		/// @param[in] timestamp    时间
		/// @return 对应的惯导参数
		const SegmentCalibParam *GetTrajParam(const double timestamp) const {
			int	segment_id = SegmentID(timestamp);
			return &calib_param_->segment_param[segment_id];
		}

		// 时间是否有效
		bool GetTrajQuality(const double timestamp) const {
			if (timestamp < this->minTime() || timestamp >= this->maxTime())
				return false;
			else
				return true;
		}

		bool GetLiDARTrajQuality(const double timestamp) const {
			double t_lidar = timestamp + this->GetTrajParam()->time_offset;

			return GetTrajQuality(t_lidar);
		}

		SE3d GetLidarPose(const double timestamp) const;

		bool GetLidarPose(const double timestamp, SE3d &lidar_pose);

		void UndistortScan(const PosCloud &scan_raw, const double target_timestamp,
			PosCloud &scan_in_target) const;

		bool EvaluateLidarRelativeRotation(double lidar_time1, double lidar_time2,
			Eigen::Quaterniond &q_L2toL1);

		void SaveTrajectoryControlPoints(std::string path);

		bool LoadTrajectoryControlPoints(std::string path);

		void TrajectoryToTUMTxt(std::string file_path, double relative_start_time = 0,
			double relative_end_time = 0, int iteration_num = 0,
			double dt = 0.02) const;

		void TrajectoryToTUMTxt2(std::string file_path,
			double relative_start_time = 0,
			double relative_end_time = 0,
			double dt = 0.02) const;

	private:
		size_t segment_id_;     ///< 轨迹分段ID
		double time_interval_;  ///< 节点的时间间隔
		CalibParamManager::Ptr calib_param_; ///< 标定参数
	};

	extern bool readGNSSPVAData(const std::string filename, std::vector<GNSSPVA> &GNSSDatas, const bool bweek = true);

	extern bool readIMUData(const std::string filename, std::vector<IMUData> &IMUDatas, const Eigen::Vector3d IMURotAngle = Eigen::Vector3d(), const int week = 0);
	extern void printOneIMUData(const IMUData &imuData, FILE *pFile = NULL);
	extern void printIMUData(const std::vector<IMUData> &IMUDatas, FILE *pFile = NULL);
	extern void syncIMUData(std::vector<IMUData> &IMUDatas, const double start_time = -1.0, const double end_time = -1.0);

	extern bool readPoseData(const int ref_type, const std::string ref_file, std::vector<PoseData> &PoseDatas, const bool bweek);
	extern void printOnePoseData(const PoseData &poseData, FILE *pFile = NULL);
	extern void printPoseData(const std::vector<PoseData> &PoseDatas, FILE *pFile = NULL);
	extern void syncPoseData(std::vector<PoseData> &PoseDatas, const double start_time = -1.0, const double end_time = -1.0);

	extern bool readOdomData(const std::string filename, std::vector<OdomData> &OdomDatas);
	extern void printOneOdomData(const OdomData &odomData, FILE *pFile = NULL);
	extern void printOdomData(const std::vector<OdomData> &OdomDatas, FILE *pFile = NULL);
	extern void syncOdomData(std::vector<OdomData> &OdomDatas, const double start_time = -1.0, const double end_time = -1.0);

	extern void PoseData2Attitude(const PoseData pose, double attitude[3], const bool bDeg);
	extern void PoseData2Attitude(const liso::SE3d pose, double attitude[3], const bool bDeg);

}  // namespace liso
#endif
