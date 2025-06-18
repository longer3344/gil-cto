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

#include <trajectory/se3_trajectory.h>
#include <algorithm>
#include "CommonFunc.h"

namespace liso {

SE3d Trajectory::GetLidarPose(const double timestamp) const {
  double t = timestamp + this->GetTrajParam()->time_offset;
  if (t < this->minTime()) {
    t = this->minTime();
  } else if (t >= this->maxTime()) {
    t = this->maxTime() - 1e-9;
  }

  SE3d pose_I_to_G = this->pose(t);
  SE3d pose_L_to_G = pose_I_to_G * calib_param_->se3_LtoI;
  return pose_L_to_G;
}

bool Trajectory::GetLidarPose(const double timestamp, SE3d &lidar_pose) {
  double t = timestamp + this->GetTrajParam(timestamp)->time_offset;
  if (t < this->minTime() || t >= this->maxTime()) return false;

  SE3d pose_I_to_G = this->pose(t);
  SE3d pose_L_to_G = pose_I_to_G * calib_param_->se3_LtoI;
  lidar_pose = pose_L_to_G;
  return true;
}

void Trajectory::UndistortScan(const PosCloud &scan_raw,
                               const double target_timestamp,
                               PosCloud &scan_in_target) const {
  scan_in_target.header = scan_raw.header;
  scan_in_target.resize(scan_raw.size());
  scan_in_target.is_dense = true;

  SE3d pose_G_to_target = GetLidarPose(target_timestamp).inverse();

  std::size_t cnt = 0;
  for (auto const &raw_p : scan_raw.points) {
    if (pcl_isnan(raw_p.x)) {
      scan_in_target.is_dense = false;
      std::cout << RED << "[UndistortScan] input cloud exists NAN point\n"
                << RESET;
      continue;
    }
    SE3d pose_Lk_to_G = GetLidarPose(raw_p.timestamp);

    Eigen::Vector3d p_Lk(raw_p.x, raw_p.y, raw_p.z);
    Eigen::Vector3d point_out;
    point_out = pose_G_to_target * pose_Lk_to_G * p_Lk;

    PosPoint point;
    point.x = point_out(0);
    point.y = point_out(1);
    point.z = point_out(2);
    point.timestamp = raw_p.timestamp;

    scan_in_target[cnt++] = point;
  }
}

bool Trajectory::EvaluateLidarRelativeRotation(double lidar_time1,
                                               double lidar_time2,
                                               Eigen::Quaterniond &q_L2toL1) {
  assert(lidar_time1 <= lidar_time2 &&
         "[EvaluateRelativeRotation] : lidar_time1 > lidar_time2");
  if (lidar_time1 < this->minTime() || lidar_time2 > this->maxTime())
    return false;

  SE3d pose_I1_to_G = this->pose(lidar_time1);
  SE3d pose_I2_to_G = this->pose(lidar_time2);
  Eigen::Quaterniond q1 = pose_I1_to_G.unit_quaternion();
  Eigen::Quaterniond q2 = pose_I2_to_G.unit_quaternion();

  Eigen::Quaterniond q_I2toI1 = q1.conjugate() * q2;
  q_L2toL1 = calib_param_->q_LtoI.conjugate() * q_I2toI1 * calib_param_->q_LtoI;
  return true;
}

void Trajectory::SaveTrajectoryControlPoints(std::string path) {
  std::ofstream outfile;
  outfile.open(path);

  size_t NumKnots = this->numKnots();
  for (size_t i = 0; i < NumKnots; i++) {
    Eigen::Vector3d p = this->getKnotPos(i);
    Sophus::SO3d s = this->getKnotSO3(i);
    Eigen::Quaterniond q = s.unit_quaternion();
    outfile << p(0) << " " << p(1) << " " << p(2) << " " << q.x() << " "
            << q.y() << " " << q.z() << " " << q.w() << "\n";
  }
  outfile.close();
}

bool Trajectory::LoadTrajectoryControlPoints(std::string path) {
  std::ifstream infile;
  infile.open(path);
  std::string current_line;

  std::vector<Eigen::VectorXd> controlPoints;
  while (std::getline(infile, current_line)) {
    std::istringstream s(current_line);
    std::string field;
    std::vector<double> vec;

    while (std::getline(s, field, ' ')) {
      if (field.empty())  // Skip if empty
        continue;
      // save the data to our vector
      vec.push_back(std::atof(field.c_str()));
    }
    // Create eigen vector
    Eigen::VectorXd temp(vec.size());
    for (size_t i = 0; i < vec.size(); i++) {
      temp(i) = vec.at(i);
    }
    controlPoints.push_back(temp);
  }

  for (unsigned int i = 0; i < controlPoints.size(); i++) {
    Eigen::VectorXd temp = controlPoints.at(i);
    Eigen::Vector3d p = temp.head<3>();
    Eigen::Quaterniond q(temp(6), temp(3), temp(4), temp(5));
    Sophus::SE3d se3_knot(q, p);
    this->knots_push_back(se3_knot);
  }

  return true;
}

void Trajectory::TrajectoryToTUMTxt(std::string file_path,
                                    double relative_start_time,
                                    double relative_end_time, int iteration_num,
                                    double dt) const {
  std::string file_name = "/trajectory-lidar-" +
                          std::to_string(relative_start_time) + "-" +
                          std::to_string(relative_end_time) + "-iter" +
                          std::to_string(iteration_num) + ".txt";
  std::string traj_path = file_path + file_name;

  TrajectoryToTUMTxt2(traj_path, relative_start_time, relative_end_time, dt);
}

void Trajectory::TrajectoryToTUMTxt2(std::string traj_path,
                                     double relative_start_time,
                                     double relative_end_time,
                                     double dt) const {
  std::ofstream outfile;
  outfile.open(traj_path);

  double min_time = this->minTime();
  double max_time = this->maxTime();
  for (double t = min_time; t < max_time; t += dt) {
    double relative_bag_time = relative_start_time + t;

    SE3d pose = this->GetLidarPose(t);
    Eigen::Vector3d p = pose.translation();
    Eigen::Quaterniond q = pose.unit_quaternion();

    outfile.precision(9);
    outfile << relative_bag_time << " ";
    outfile.precision(5);
    outfile << p(0) << " " << p(1) << " " << p(2) << " " << q.x() << " "
            << q.y() << " " << q.z() << " " << q.w() << "\n";
  }
  outfile.close();
  std::cout << "Save trajectory at " << traj_path << std::endl;
}

bool readGNSSPVAData(const std::string filename, std::vector<GNSSPVA> &GNSSDatas, const bool bweek)
{
	char    buff[1024] = { '\0' };
	FILE	*POS_fp = fopen(filename.c_str(), "r");
	if (POS_fp == NULL)
	{
		GNSSDatas.clear();
		return false;
	}

	int week = 0;
	double sec = 0.0;
	ELC_PVA     POS_DATA;
	GNSSPVA		gnss;

	GNSSDatas.clear();
	while (!feof(POS_fp))
	{
		if (fgets(buff, sizeof(buff), POS_fp) == NULL)
			break;

		memset(&POS_DATA, 0, sizeof(POS_DATA));
		int num = sscanf(buff, "%d %lf %hhd %hhd %f %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			&(week), &(sec),
			&(POS_DATA.Q), &(POS_DATA.AmbFix), &(POS_DATA.DDOP),
			POS_DATA.Pos + 0, POS_DATA.Pos + 1, POS_DATA.Pos + 2,
			POS_DATA.PosP + 0, POS_DATA.PosP + 1, POS_DATA.PosP + 2,
			POS_DATA.PosP + 3, POS_DATA.PosP + 4, POS_DATA.PosP + 5,
			POS_DATA.Vel + 0, POS_DATA.Vel + 1, POS_DATA.Vel + 2,
			POS_DATA.VelP + 0, POS_DATA.VelP + 1, POS_DATA.VelP + 2,
			POS_DATA.VelP + 3, POS_DATA.VelP + 4, POS_DATA.VelP + 5,
			POS_DATA.Att + 0, POS_DATA.Att + 1, POS_DATA.Att + 2,
			POS_DATA.AttBL + 0, POS_DATA.AttBL + 1, POS_DATA.AttBL + 2,
			POS_DATA.AttP + 0, POS_DATA.AttP + 1, POS_DATA.AttP + 2);

		if (num < 11) continue;

		gnss.week = week;
		gnss.sow = sec;
		gnss.timestamp = bweek ? week * 7 * 86400.0 + sec : sec;

		for (int i = 0; i < 3; i++)
		{
			gnss.position(i) = POS_DATA.Pos[i];
			gnss.velocity(i) = POS_DATA.Vel[i];
		}

		gnss.poscov.Zero();
		gnss.poscov(0, 0) = POS_DATA.PosP[0];  gnss.poscov(0, 1) = POS_DATA.PosP[3];  gnss.poscov(0, 2) = POS_DATA.PosP[4];
		gnss.poscov(1, 0) = POS_DATA.PosP[3];  gnss.poscov(1, 1) = POS_DATA.PosP[1];  gnss.poscov(1, 2) = POS_DATA.PosP[5];
		gnss.poscov(2, 0) = POS_DATA.PosP[4];  gnss.poscov(2, 1) = POS_DATA.PosP[5];  gnss.poscov(2, 2) = POS_DATA.PosP[2];

		gnss.velcov.Zero();
		gnss.velcov(0, 0) = POS_DATA.VelP[0];  gnss.velcov(0, 1) = POS_DATA.VelP[3];  gnss.velcov(0, 2) = POS_DATA.VelP[4];
		gnss.velcov(1, 0) = POS_DATA.VelP[3];  gnss.velcov(1, 1) = POS_DATA.VelP[1];  gnss.velcov(1, 2) = POS_DATA.VelP[5];
		gnss.velcov(2, 0) = POS_DATA.VelP[4];  gnss.velcov(2, 1) = POS_DATA.VelP[5];  gnss.velcov(2, 2) = POS_DATA.VelP[2];


		//std::cout << gnss.timestamp << std::endl;
		//std::cout << gnss.poscov << std::endl;
		GNSSDatas.push_back(gnss);
	}
	fclose(POS_fp);


	for (size_t i = 0; i < GNSSDatas.size(); i++)
	{
		if (std::isinf(GNSSDatas[i].position(0)) || std::isinf(GNSSDatas[i].position(1)) || std::isinf(GNSSDatas[i].position(2)))
		{
			int xxx = 0;
		}

		if (std::isinf(GNSSDatas[i].poscov(0, 0)) || std::isinf(GNSSDatas[i].poscov(1, 1)) || std::isinf(GNSSDatas[i].poscov(2, 2)))
		{
			int xxx = 0;
		}
	}

	if (GNSSDatas.size() > 0)
		printf("GNSS Data : %10zd %10.3lf %10.3lf\n", GNSSDatas.size(), GNSSDatas.front().sow, GNSSDatas.back().sow);
	return GNSSDatas.size() > 0;
}

bool readIMUData(const std::string filename, std::vector<IMUData> &IMUDatas, const Eigen::Vector3d IMURotAngle, const int week)
{
	IMUDatas.clear();
	Eigen::Matrix3d		Rpb = Eigen::Matrix3d(Eigen::AngleAxisd(-IMURotAngle[1], Eigen::Vector3d::UnitY()) *
		Eigen::AngleAxisd(-IMURotAngle[0], Eigen::Vector3d::UnitX()) *
		Eigen::AngleAxisd(-IMURotAngle[2], Eigen::Vector3d::UnitZ()));
	Eigen::Matrix3d		Rbp = Rpb.transpose();
	std::cout << "Rbp" << std::endl;
	std::cout << Rbp << std::endl;

	FILE*	fp = NULL;
	fp = fopen(filename.c_str(), "r");
	if (!fp)	return false;

	IMUData	imu;
	char	oneline[256] = { 0 };
	Eigen::Vector3d gyro;
	Eigen::Vector3d accel;
	double	sow = 0.0;
	while (!feof(fp))
	{
		if (fgets(oneline, sizeof(oneline), fp) == NULL)
			break;

		if (sscanf(oneline, "%lf %lf %lf %lf %lf %lf %lf", &sow, gyro.data() + 0, gyro.data() + 1, gyro.data() + 2, 
			accel.data() + 0, accel.data() + 1, accel.data() + 2) != 7)
			continue;

		gyro = gyro * D2R;
		imu.timestamp = week * 7 * 86400.0 + sow;
		imu.gyro = Rbp * gyro;
		imu.accel = Rbp * accel;
		imu.orientation.setIdentity();

		//printOneIMUData(imu);
		IMUDatas.push_back(imu);
	}
	fclose(fp);
	
	if (IMUDatas.size() > 0)
		printf("IMU  Data : %10zd %10.3lf %10.3lf\n", IMUDatas.size(), IMUDatas.front().timestamp, IMUDatas.back().timestamp);
	return IMUDatas.size() > 0;
}
void printOneIMUData(const IMUData &imuData, FILE *pFile)
{
	FILE	*Fout = pFile ? pFile : stdout;
	double	attitude[3] = { 0.0 }
	;
	Eigen::Quaterniond qbn = imuData.orientation;
	Eigen::Matrix3d Rbn = qbn.toRotationMatrix();
	RotMatrix2RotAngle(Rbn.transpose().data(), attitude);

	fprintf(Fout, "%12.6lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %9.6lf %9.6lf %9.6lf %9.6lf %11.6lf %10.6lf %10.6lf\n", fmod(imuData.timestamp, 604800.0),
		imuData.gyro(0) * R2D, imuData.gyro(1) * R2D, imuData.gyro(2) * R2D,
		imuData.accel(0), imuData.accel(1), imuData.accel(2),
		qbn.x(), qbn.y(), qbn.z(), qbn.w(), attitude[0] * R2D, attitude[1] * R2D, attitude[2] * R2D);
}
void printIMUData(const std::vector<IMUData> &IMUDatas, FILE *pFile)
{
	for (size_t i = 0; i < IMUDatas.size(); i++)
		printOneIMUData(IMUDatas[i]);
	return;
}
void syncIMUData(std::vector<IMUData> &IMUDatas, const double start_time, const double end_time)
{
	if (IMUDatas.size() <= 0 || (start_time < 0 && end_time < 0))
		return;

	if (start_time >= 0)
	{
		while (IMUDatas.size() > 0)
		{
			if (fabs(start_time - IMUDatas.front().timestamp) < 1.0e-06)
				break;

			IMUDatas.erase(std::begin(IMUDatas));
		}
	}

	if (end_time >= 0)
	{
		while (IMUDatas.size() > 0)
		{
			if (fabs(end_time - IMUDatas.back().timestamp) < 1.0e-06)
				break;

			IMUDatas.pop_back();
		}
	}

	return;
}

bool readPoseData(const int ref_type, const std::string ref_file, std::vector<PoseData> &PoseDatas, const bool bweek)
{
	PoseDatas.clear();

	FILE*	fp = NULL;
	fp = fopen(ref_file.c_str(), "r");
	if (!fp)	return false;

	PoseData	pose;
	char	oneline[1024] = { 0 };
	Eigen::Vector3d gyro;
	Eigen::Vector3d accel;
	int		week = 0;
	double	sow = 0.0, attitude[3] = { 0.0 }, LLH[3] = { 0.0 }, Rnb[9] = { 0.0 }, Rbe[9] = { 0.0 }, Reb[9] = { 0.0 };
	while (!feof(fp))
	{
		if (fgets(oneline, sizeof(oneline), fp) == NULL)
			break;

		if (ref_type == 0)
		{
			// IPS
			if (sscanf(oneline, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &week, &sow, 
				pose.position.data() + 0, pose.position.data() + 1, pose.position.data() + 2,
				pose.velocity.data() + 0, pose.velocity.data() + 1, pose.velocity.data() + 2,
				attitude + 0, attitude + 1, attitude + 2) != 11)
				continue;
		}
		else if (ref_type == 1)
		{
			// IE
			if (sscanf(oneline, "%d %lf %*f %*f %*f %*d %*s %*s %*s %lf %lf %lf %*f %*f %*f %lf %lf %lf %*f %*f %*f %lf %lf %lf", &week, &sow,
				pose.position.data() + 0, pose.position.data() + 1, pose.position.data() + 2,
				pose.velocity.data() + 0, pose.velocity.data() + 1, pose.velocity.data() + 2,
				attitude + 0, attitude + 1, attitude + 2) != 11)
				continue;
		}
		else	continue;

		pose.week = week; pose.sow = sow;
		pose.timestamp = bweek ? week * 7 * 86400.0 + sow : sow;
		XYZ2LLH(pose.position, LLH);
		attitude[0] *= D2R; attitude[1] *= D2R; attitude[2] *= D2R;
		Azimuth2Attitude(attitude, attitude);
		Attitude2Rbe(attitude, LLH, Rbe);
		MatrixTranspose(3, 3, Rbe, Reb);
		RotAngle2RotMatrix(attitude, Rnb);

		Eigen::Matrix3d rot;
		rot << Rbe[0], Rbe[1], Rbe[2], Rbe[3], Rbe[4], Rbe[5], Rbe[6], Rbe[7], Rbe[8];
		//std::cout << rot << std::endl << std::endl;

		Eigen::Quaterniond quat = Eigen::Quaterniond(rot);
		pose.orientation = liso::SO3d(quat);

		//printOnePoseData(pose);
		PoseDatas.push_back(pose);
	}
	fclose(fp);
	
	if (PoseDatas.size() > 0)
		printf("Pose Data : %10zd %10.3lf %10.3lf\n", PoseDatas.size(), PoseDatas.front().sow, PoseDatas.back().sow);
	return PoseDatas.size() > 0;
}

void printOnePoseData(const PoseData &poseData, FILE *pFile)
{
	FILE	*Fout = pFile ? pFile : stdout;
	double	LLH[3] = { 0.0 }, Rbe[9] = { 0.0 }, attitude[3] = { 0.0 }, azimuth[3] = { 0.0 };

	Eigen::Quaterniond qbn = poseData.orientation.unit_quaternion();
	Eigen::Matrix3d Rbn = qbn.toRotationMatrix();
	XYZ2LLH(poseData.position.data(), LLH);
	//liso::Rbe2Attitude(Rbn.data(), LLH, attitude);
	Rbe[0] = Rbn(0, 0); Rbe[1] = Rbn(0, 1); Rbe[2] = Rbn(0, 2);
	Rbe[3] = Rbn(1, 0); Rbe[4] = Rbn(1, 1); Rbe[5] = Rbn(1, 2);
	Rbe[6] = Rbn(2, 0); Rbe[7] = Rbn(2, 1); Rbe[8] = Rbn(2, 2);
	Rbe2Attitude(Rbe, LLH, attitude);
	Attitude2Azimuth(attitude, azimuth);

	fprintf(Fout, "%4d %12.5lf %14.5lf %14.5lf %14.5lf %9.6lf %9.6lf %9.6lf %9.6lf %9.6lf %9.6lf %9.6lf %11.6lf %10.6lf %10.6lf\n", int(poseData.timestamp / 604800.0), fmod(poseData.timestamp, 604800.0),
		poseData.position(0), poseData.position(1), poseData.position(2),
		poseData.velocity(0), poseData.velocity(1), poseData.velocity(2),
		qbn.x(), qbn.y(), qbn.z(), qbn.w(), azimuth[0] * R2D, azimuth[1] * R2D, azimuth[2] * R2D);

	return;
}
void printPoseData(const std::vector<PoseData> &PoseDatas, FILE *pFile)
{
	for (size_t i = 0; i < PoseDatas.size(); i++)
		printOnePoseData(PoseDatas[i]);
	return;
}
void syncPoseData(std::vector<PoseData> &PoseDatas, const double start_time, const double end_time)
{
	if (PoseDatas.size() <= 0 || (start_time < 0 && end_time < 0))
		return;

	if (start_time >= 0)
	{
		while (PoseDatas.size() > 0)
		{
			if (fabs(start_time - PoseDatas.front().timestamp) < 1.0e-06)
				break;

			PoseDatas.erase(std::begin(PoseDatas));
		}
	}

	if (end_time >= 0)
	{
		while (PoseDatas.size() > 0)
		{
			if (fabs(end_time - PoseDatas.back().timestamp) < 1.0e-06)
				break;

			PoseDatas.pop_back();
		}
	}

	return;
}

bool readOdomData(const std::string filename, std::vector<OdomData> &OdomDatas)
{
	OdomDatas.clear();

	FILE*	fp = NULL;
	fp = fopen(filename.c_str(), "r");
	if (!fp)	return false;


	fclose(fp);
	return OdomDatas.size() > 0;
}
void printOneOdomData(const OdomData &odomData, FILE *pFile)
{
	FILE	*Fout = pFile ? pFile : stdout;
	double	LLH[3] = { 0.0 }, attitude[3] = { 0.0 }, azimuth[3] = { 0.0 };

	Eigen::Matrix3d Rbn = odomData.pose.block(0, 0, 3, 3);
	Eigen::Quaterniond qbn = Eigen::Quaterniond(Rbn);
	RotMatrix2RotAngle(Rbn.transpose().data(), attitude);

	fprintf(Fout, "%12.5lf %14.5lf %14.5lf %14.5lf %9.6lf %9.6lf %9.6lf %9.6lf %9.6lf %9.6lf %9.6lf\n", odomData.timestamp,
		odomData.pose(0, 3), odomData.pose(1, 3), odomData.pose(2, 3),
		qbn.x(), qbn.y(), qbn.z(), qbn.w(), azimuth[0], azimuth[1], azimuth[2]);

	return;
}
void printOdomData(const std::vector<OdomData> &OdomDatas, FILE *pFile)
{
	for (size_t i = 0; i < OdomDatas.size(); i++)
		printOneOdomData(OdomDatas[i]);
	return ;
}
void syncOdomData(std::vector<OdomData> &OdomDatas, const double start_time, const double end_time)
{
	if (OdomDatas.size() <= 0 || (start_time < 0 && end_time < 0))
		return;

	if (start_time >= 0)
	{
		while (OdomDatas.size() > 0)
		{
			if (fabs(start_time - OdomDatas.front().timestamp) < 1.0e-06)
				break;

			OdomDatas.erase(std::begin(OdomDatas));
		}
	}

	if (end_time >= 0)
	{
		while (OdomDatas.size() > 0)
		{
			if (fabs(end_time - OdomDatas.back().timestamp) < 1.0e-06)
				break;

			OdomDatas.pop_back();
		}
	}

	return;
}


void PoseData2Attitude(const PoseData pose, double attitude[3], const bool bDeg)
{
	double	Rbe[9], tbe[3], LLH[3];

	Eigen::Matrix3d rr = pose.orientation.unit_quaternion().matrix();
	Rbe[0] = rr(0, 0); Rbe[1] = rr(0, 1); Rbe[2] = rr(0, 2);
	Rbe[3] = rr(1, 0); Rbe[4] = rr(1, 1); Rbe[5] = rr(1, 2);
	Rbe[6] = rr(2, 0); Rbe[7] = rr(2, 1); Rbe[8] = rr(2, 2);
	tbe[0] = pose.position[0]; tbe[1] = pose.position[1]; tbe[2] = pose.position[2];
	XYZ2LLH(pose.position, LLH);
	Rbe2Attitude(Rbe, LLH, attitude);
	Attitude2Azimuth(attitude, attitude);
	if (bDeg)
	{
		attitude[0] *= R2D; attitude[1] *= R2D; attitude[2] *= R2D;
	}

	return;
}

void PoseData2Attitude(const liso::SE3d pose, double attitude[3], const bool bDeg)
{
	double	Rbe[9], tbe[3], LLH[3];

	Eigen::Matrix3d rr = pose.unit_quaternion().toRotationMatrix();
	Rbe[0] = rr(0, 0); Rbe[1] = rr(0, 1); Rbe[2] = rr(0, 2);
	Rbe[3] = rr(1, 0); Rbe[4] = rr(1, 1); Rbe[5] = rr(1, 2);
	Rbe[6] = rr(2, 0); Rbe[7] = rr(2, 1); Rbe[8] = rr(2, 2);
	tbe[0] = pose.translation()[0]; tbe[1] = pose.translation()[1]; tbe[2] = pose.translation()[2];
	XYZ2LLH(tbe, LLH);
	Rbe2Attitude(Rbe, LLH, attitude);
	Attitude2Azimuth(attitude, attitude);
	if (bDeg)
	{
		attitude[0] *= R2D; attitude[1] *= R2D; attitude[2] *= R2D;
	}

	return;
}

}  // namespace liso
