#include <stdio.h>

#include "basalt/spline/spline_segment.h"
#include "trajectory/se3_trajectory.h"
#include "trajectory/trajectory_estimator.h"

#include "LiDARFFStream.h"

#include <ceres/ceres.h>

#include "HC_Utils.hpp"

#include <basalt/spline/rd_spline.h>
#include <basalt/spline/so3_spline.h>

static double ValidDouble(const double value, const double eps)
{
	double	scale = 1.0 / eps;
	long long  lval = (long long)(value * scale + 0.5);
	double ret = lval * eps;
	return ret;
}


void SE3dToAzimuth(liso::SE3d pose, double Azimuth[3])
{
	double	LLH[3] = { 0.0 }, Rbe[9] = { 0.0 }, attitude[3] = { 0.0 };
	
	XYZ2LLH(pose.translation(), LLH);
	Eigen::Matrix3d rr = pose.unit_quaternion().toRotationMatrix();
	Rbe[0] = rr(0, 0); Rbe[1] = rr(0, 1); Rbe[2] = rr(0, 2);
	Rbe[3] = rr(1, 0); Rbe[4] = rr(1, 1); Rbe[5] = rr(1, 2);
	Rbe[6] = rr(2, 0); Rbe[7] = rr(2, 1); Rbe[8] = rr(2, 2);

	Rbe2Attitude(Rbe, LLH, Azimuth);
	Attitude2Azimuth(Azimuth, Azimuth);
	Azimuth[0] *= R2D; Azimuth[1] *= R2D; Azimuth[2] *= R2D;
}

void PoseDataToAzimuth(liso::PoseData pose, double Azimuth[3])
{
	double	LLH[3] = { 0.0 }, Rbe[9] = { 0.0 }, attitude[3] = { 0.0 };

	XYZ2LLH(pose.position, LLH);
	Eigen::Matrix3d rr = pose.orientation.unit_quaternion().toRotationMatrix();
	Rbe[0] = rr(0, 0); Rbe[1] = rr(0, 1); Rbe[2] = rr(0, 2);
	Rbe[3] = rr(1, 0); Rbe[4] = rr(1, 1); Rbe[5] = rr(1, 2);
	Rbe[6] = rr(2, 0); Rbe[7] = rr(2, 1); Rbe[8] = rr(2, 2);

	Rbe2Attitude(Rbe, LLH, Azimuth);
	Attitude2Azimuth(Azimuth, Azimuth);
	Azimuth[0] *= R2D; Azimuth[1] *= R2D; Azimuth[2] *= R2D;
}


template <int N>
void test_ceres_spline_helper_so3() {
	static const double dt = 2.0;

	basalt::So3Spline<N> spline(dt);
	spline.genRandomTrajectory(3 * N);

	for (double t = 0; t < spline.maxTime(); t += 0.1) {
		Sophus::SO3d pos1 = spline.evaluate(t);
		Eigen::Vector3d vel1 = spline.velocityBody(t);
		Eigen::Vector3d accel1 = spline.accelerationBody(t);
		Eigen::Vector3d jerk1 = spline.jerkBody(t);

		Eigen::Quaterniond q1 = pos1.unit_quaternion();
		std::cout << "t = " << t << std::endl;
		std::cout << "q1     : " << q1.x() << q1.y() << q1.z() << q1.w() << std::endl;
		std::cout << "vel1   : " << vel1.transpose() << std::endl;
		std::cout << "accel1 : " << accel1.transpose() << std::endl;
		std::cout << "jerk1  : " << jerk1.transpose() << std::endl;

		Sophus::SO3d pos2;
		Eigen::Vector3d vel2, accel2, jerk2;

		{
			double pow_inv_dt = 1.0 / dt;

			int64_t st_ns = (t);

			BASALT_ASSERT_STREAM(st_ns >= 0, "st_ns " << st_ns << " time_ns " << t_ns
				<< " start_t_ns " << 0);

			int64_t s = std::floor(t / dt);
			double u = t - s*dt;

			BASALT_ASSERT_STREAM(s >= 0, "s " << s);
			BASALT_ASSERT_STREAM(size_t(s + N) <= spline.getKnots().size(),
				"s " << s << " N " << N << " knots.size() "
				<< spline.getKnots().size());

			std::vector<const double *> vec;
			for (int i = 0; i < N; i++) {
				vec.emplace_back(spline.getKnots()[s + i].data());
			}

			basalt::CeresSplineHelper<N>::template evaluate_lie<double, Sophus::SO3>(
				&vec[0], u, pow_inv_dt, &pos2);
			basalt::CeresSplineHelper<N>::template evaluate_lie<double, Sophus::SO3>(
				&vec[0], u, pow_inv_dt, nullptr, &vel2);
			basalt::CeresSplineHelper<N>::template evaluate_lie<double, Sophus::SO3>(
				&vec[0], u, pow_inv_dt, nullptr, nullptr, &accel2);
			basalt::CeresSplineHelper<N>::template evaluate_lie<double, Sophus::SO3>(
				&vec[0], u, pow_inv_dt, nullptr, nullptr, nullptr, &jerk2);

			Eigen::Quaterniond q2 = pos1.unit_quaternion();
			std::cout << "q2     : " << q2.x() << q2.y() << q2.z() << q2.w() << std::endl;
			std::cout << "vel2   : " << vel2.transpose() << std::endl;
			std::cout << "accel2 : " << accel2.transpose() << std::endl;
			std::cout << "jerk2  : " << jerk2.transpose() << std::endl;
		}
	}
}


long unsigned int MapPoint::nNextId = 0;	///< 点的全局ID
long unsigned int MapLine::nNextId = 0;		///< 线的全局ID
long unsigned int MapPlane::nNextId = 0;	///< 面的全局ID


static bool initTrajectory(std::vector<liso::PoseData> &posedata, std::shared_ptr<liso::Trajectory> &trajectory)
{
	// 权重参数
	double	pos_weight = 5.0;	// 位置权重
	double	vel_weight = 1.0;	// 速度权重
	double	rot_weight = 10.0;	// 姿态权重
	double	gyr_weight = 1.0 / 0.118; // 陀螺权重
	double	acc_weight = 1.0 / 0.664; // 加计权重

	liso::TrajectoryEstimatorOptions options;
	options.lock_ab = true;	    // 是否不优化加计零偏
	options.lock_wb = true;	    // 是否不优化陀螺零偏
	options.lock_g = true;	    // 是否不优化重力
	options.lock_tlb = true;    // 是否不优化杆臂
	options.lock_Rlb = true;    // 是否不优化安置角
	options.lock_trajP = false; // 是否不优化位置节点
	options.lock_trajR = false; // 是否不优化姿态节点
	options.local_traj = false; // 非局部系，直接用位置计算重力

	liso::TrajectoryEstimator estimator(trajectory, nullptr, options);
	for (int i = 0; i < posedata.size(); i++)
	{
		//printf("%04d %.8lf %.8lf\n", posedata[i].week, posedata[i].sow, posedata[i].timestamp);
		estimator.AddPositionMeasurement(posedata[i], pos_weight);		// 位置
		estimator.AddVelocityMeasurement(posedata[i], vel_weight);		// 速度
		estimator.AddOrientationMeasurement(posedata[i], rot_weight);	// 姿态
	}

	printf("initTrajectory ParameterBlocks = %d\n", estimator.problem_->NumParameterBlocks());
	printf("initTrajectory ResidualBlocks  = %d\n", estimator.problem_->NumResidualBlocks());
	ceres::Solver::Summary summary = estimator.Solve(50, true);
	std::cout << "[initTrajectory]: " << summary.BriefReport() << std::endl;

	FILE *fout = fopen("initTrajectory.txt", "w");
	for (size_t i = 0; i < posedata.size(); i++)
	{
		liso::PoseData pose0 = posedata[i];
		liso::SE3d pose = trajectory->pose(posedata[i].timestamp);
		Eigen::Vector3d dpos = pose.translation() - posedata[i].position;

		if (fabs(posedata[i].timestamp - 553938.000) < 0.001)
			printf("");

		double	Azimuth0[3] = { 0.0 }, Azimuth1[3] = { 0.0 };
		Sophus::SO3d::Tangent dPhi = (pose.so3()*pose0.orientation.inverse()).log();
		PoseDataToAzimuth(posedata[i], Azimuth0);
		SE3dToAzimuth(pose, Azimuth1);
		fprintf(fout, "%13.6lf %8.3lf %8.3lf %8.3lf %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf %11.6lf\n", posedata[i].timestamp,
			dpos(0), dpos(1), dpos(2), Azimuth1[0] - Azimuth0[0], Azimuth1[1] - Azimuth0[1], Azimuth1[2] - Azimuth0[2],
			dPhi[0], dPhi[1], dPhi[2]);
	}
	if (fout) { fclose(fout); fout = NULL; }

	return true;
}

static bool initIMUAccelBias(std::vector<liso::IMUData> &IMUdata,
	std::shared_ptr<liso::Trajectory> &trajectory, liso::CalibParamManager::Ptr &calib_param)
{
	double	pos_weight = 1.0;
	double	vel_weight = 5.0;
	double	rot_weight = 10.0;
	double	gyr_weight = 1.0 / 0.118;
	double	acc_weight = 1.0 / 0.664;
	double	GyrPNSD = 5.06797954e-06;
	double	AccPNSD = 9.00000000e-05;

	liso::TrajectoryEstimatorOptions options;
	options.lock_ab = false;    // 是否不优化加计零偏
	options.lock_wb = false;    // 是否不优化陀螺零偏
	options.lock_g = true;      // 是否不优化重力
	options.lock_tlb = true;    // 是否不优化杆臂
	options.lock_Rlb = true;    // 是否不优化安置角
	options.lock_trajP = false; // 是否不优化位置节点
	options.lock_trajR = false; // 是否不优化姿态节点
	options.local_traj = false; // 非局部系，直接用位置计算重力

	liso::TrajectoryEstimator estimator(trajectory, calib_param, options);
	for (int i = 0; i < IMUdata.size(); i++)
	{
		//printf("%04d %.8lf %.8lf\n", IMUdata[i].timestamp, IMUdata[i].accel(0), IMUdata[i].accel(1), IMUdata[i].accel(2));
		estimator.AddIMUAccelMeasurement(IMUdata[i], acc_weight);

		if (fabs(std::remainder(IMUdata[i].timestamp, 1.0)) < 0.003)
		{
			Eigen::Vector3d acc_residual = trajectory->accelResidual(IMUdata[i].timestamp, IMUdata[i].accel);
			Eigen::Vector3d gyr_residual = trajectory->gyroResidual(IMUdata[i].timestamp, IMUdata[i].gyro);
			fprintf(stdout, "%10.3lf %14.5lf %14.5lf %14.5lf %14.5lf %14.5lf %14.5lf\n", IMUdata[i].timestamp,
				acc_residual(0), acc_residual(1), acc_residual(2), gyr_residual(0), gyr_residual(1), gyr_residual(2));
		}
	}

	// 先验约束
	//estimator.AddAccelPriorBias(5E-02);

	// 随机游走约束
	//estimator.AddAccelBiasBetweenConstraint(AccPNSD);

	printf("initIMUAccelBias ParameterBlocks = %d\n", estimator.problem_->NumParameterBlocks());
	printf("initIMUAccelBias ResidualBlocks  = %d\n", estimator.problem_->NumResidualBlocks());
	ceres::Solver::Summary summary = estimator.Solve(50, true);
	std::cout << "[initIMUAccelBias]: " << summary.BriefReport() << std::endl;

	for (size_t i = 0; i < IMUdata.size(); i++)
	{
		//printf("%04d %.8lf %.8lf\n", IMUdata[i].timestamp, IMUdata[i].accel(0), IMUdata[i].accel(1), IMUdata[i].accel(2));
		auto param = trajectory->GetTrajParam(IMUdata[i].timestamp); // 根据时间选取对应时间段的重力和零偏
		if (fabs(std::remainder(IMUdata[i].timestamp, 1.0)) < 0.003)
		{
			Eigen::Vector3d acc_residual = trajectory->accelResidual(IMUdata[i].timestamp, IMUdata[i].accel);
			Eigen::Vector3d gyr_residual = trajectory->gyroResidual(IMUdata[i].timestamp, IMUdata[i].gyro);
			fprintf(stdout, "%10.3lf %14.5lf %14.5lf %14.5lf %14.5lf %14.5lf %14.5lf\n", IMUdata[i].timestamp,
				acc_residual(0), acc_residual(1), acc_residual(2), gyr_residual(0), gyr_residual(1), gyr_residual(2));
		}
	}

	calib_param->showIMUBias();



	//FILE *fout = fopen("initIMUAccelBias.txt", "w");
	//for (size_t i = 0; i < calib_param->segment_param.size(); i++)
	//{
	//	fprintf(fout, "%8zd %10.3lf %10.3lf %13.10lf %13.10lf %13.10lf %13.10lf %13.10lf %13.10lf\n", i,
	//		calib_param->segment_param[i].starttime, calib_param->segment_param[i].endtime, 
	//		calib_param->segment_param[i].acce_bias(0), calib_param->segment_param[i].acce_bias(1), calib_param->segment_param[i].acce_bias(2),
	//		calib_param->segment_param[i].gyro_bias(0), calib_param->segment_param[i].gyro_bias(1), calib_param->segment_param[i].gyro_bias(2));
	//}
	//if (fout) { fclose(fout); fout = NULL; }

	return true;
}

static bool initIMUGyrosBias(std::vector<liso::IMUData> &IMUdata,
	std::shared_ptr<liso::Trajectory> &trajectory, liso::CalibParamManager::Ptr &calib_param)
{
	double	pos_weight = 1.0;
	double	vel_weight = 5.0;
	double	rot_weight = 10.0;
	double	gyr_weight = 1.0 / 0.118;
	double	acc_weight = 1.0 / 0.664;
	double	GyrPNSD = 8.00000000e-05;
	double	AccPNSD = 9.00000000e-05;

	liso::TrajectoryEstimatorOptions options;
	options.lock_ab = false;    // 是否不优化加计零偏
	options.lock_wb = false;    // 是否不优化陀螺零偏
	options.lock_g = true;      // 是否不优化重力
	options.lock_tlb = true;    // 是否不优化杆臂
	options.lock_Rlb = true;    // 是否不优化安置角
	options.lock_trajP = false; // 是否不优化位置节点
	options.lock_trajR = false; // 是否不优化姿态节点
	options.local_traj = false; // 非局部系，直接用位置计算重力

	liso::TrajectoryEstimator estimator(trajectory, calib_param, options);
	for (int i = 0; i < IMUdata.size(); i++)
	{
		//printf("%04d %.8lf %.8lf\n", IMUdata[i].timestamp, IMUdata[i].accel(0), IMUdata[i].accel(1), IMUdata[i].accel(2));
		estimator.AddIMUGyroMeasurement(IMUdata[i], gyr_weight);

		if (fabs(std::remainder(IMUdata[i].timestamp, 1.0)) < 0.003)
		{
			Eigen::Vector3d acc_residual = trajectory->accelResidual(IMUdata[i].timestamp, IMUdata[i].accel);
			Eigen::Vector3d gyr_residual = trajectory->gyroResidual(IMUdata[i].timestamp, IMUdata[i].gyro);
			fprintf(stdout, "%10.3lf %14.5lf %14.5lf %14.5lf %14.5lf %14.5lf %14.5lf\n", IMUdata[i].timestamp,
				acc_residual(0), acc_residual(1), acc_residual(2), gyr_residual(0), gyr_residual(1), gyr_residual(2));
		}
	}

	// 先验约束
	//estimator.AddGyroPriorBias(250*D2R/3600.0);

	// 随机游走约束
	//estimator.AddGyroBiasBetweenConstraint(GyrPNSD);

	printf("ParameterBlocks = %d\n", estimator.problem_->NumParameterBlocks());
	printf("ResidualBlocks  = %d\n", estimator.problem_->NumResidualBlocks());
	ceres::Solver::Summary summary = estimator.Solve(50, true);
	std::cout << "[initTrajectory]: " << summary.BriefReport() << std::endl;

	for (size_t i = 0; i < IMUdata.size(); i++)
	{
		//printf("%04d %.8lf %.8lf\n", IMUdata[i].timestamp, IMUdata[i].accel(0), IMUdata[i].accel(1), IMUdata[i].accel(2));
		auto param = trajectory->GetTrajParam(IMUdata[i].timestamp); // 根据时间选取对应时间段的重力和零偏
		if (fabs(std::remainder(IMUdata[i].timestamp, 1.0)) < 0.003)
		{
			Eigen::Vector3d acc_residual = trajectory->accelResidual(IMUdata[i].timestamp, IMUdata[i].accel, param->acce_bias);
			Eigen::Vector3d gyr_residual = trajectory->gyroResidual(IMUdata[i].timestamp, IMUdata[i].gyro, param->gyro_bias);
			fprintf(stdout, "%10.3lf %14.5lf %14.5lf %14.5lf %14.5lf %14.5lf %14.5lf\n", IMUdata[i].timestamp,
				acc_residual(0), acc_residual(1), acc_residual(2), gyr_residual(0), gyr_residual(1), gyr_residual(2));
		}
	}

	calib_param->showIMUBias();

	return true;
}

static bool initIMUBias(std::vector<liso::IMUData> &IMUdata, 
	std::shared_ptr<liso::Trajectory> &trajectory, liso::CalibParamManager::Ptr &calib_param)
{
	double	pos_weight = 1.0;
	double	vel_weight = 5.0;
	double	rot_weight = 10.0;
	double	gyr_weight = 1.0 / 0.118;
	double	acc_weight = 1.0 / 0.664;

	liso::TrajectoryEstimatorOptions options;
	options.lock_ab = false;    // 是否不优化加计零偏
	options.lock_wb = false;    // 是否不优化陀螺零偏
	options.lock_g = true;      // 是否不优化重力
	options.lock_tlb = true;    // 是否不优化杆臂
	options.lock_Rlb = true;    // 是否不优化安置角
	options.lock_trajP = false; // 是否不优化位置节点
	options.lock_trajR = false; // 是否不优化姿态节点
	options.local_traj = false; // 非局部系，直接用位置计算重力

	liso::TrajectoryEstimator estimator(trajectory, calib_param, options);
	for (int i = 0; i < IMUdata.size(); i++)
	{
		//printf("%04d %.8lf %.8lf\n", IMUdata[i].timestamp, IMUdata[i].accel(0), IMUdata[i].accel(1), IMUdata[i].accel(2));
		estimator.AddIMUAccelMeasurement(IMUdata[i], acc_weight);
		estimator.AddIMUGyroMeasurement(IMUdata[i], gyr_weight);
	}

	printf("ParameterBlocks = %d\n", estimator.problem_->NumParameterBlocks());
	printf("ResidualBlocks  = %d\n", estimator.problem_->NumResidualBlocks());
	ceres::Solver::Summary summary = estimator.Solve(50, true);
	std::cout << "[initTrajectory]: " << summary.BriefReport() << std::endl;

	for (size_t i = 0; i < IMUdata.size(); i++)
	{
		//printf("%04d %.8lf %.8lf\n", IMUdata[i].timestamp, IMUdata[i].accel(0), IMUdata[i].accel(1), IMUdata[i].accel(2));
		auto param = trajectory->GetTrajParam(IMUdata[i].timestamp); // 根据时间选取对应时间段的重力和零偏
		Eigen::Vector3d acc_residual = trajectory->accelResidual(IMUdata[i].timestamp, IMUdata[i].accel, param->acce_bias);
		Eigen::Vector3d gyr_residual = trajectory->gyroResidual(IMUdata[i].timestamp, IMUdata[i].gyro, param->gyro_bias);
	}

	return true;
}


static void get_imu_residual()
{
	// 读取位姿
	int		ref_type = 1;
	std::string	ref_file = "C:\\Users\\Junlong\\Desktop\\2023-01-05-01\\363880-363990\\ROVE_20230105_01_STIM-TC-CS.imr.ref";
	std::vector<liso::PoseData> posedata;
	liso::readPoseData(ref_type, ref_file, posedata, false);

	const double traj_start = posedata[0].timestamp;
	const double pose_interval = 0.1;
	const double traj_end = posedata.back().timestamp;

	// 初始化轨迹 knot_distance - 分段时间间隔
	std::shared_ptr<liso::Trajectory> trajectory = std::make_shared<liso::Trajectory>(pose_interval, traj_start, 0);
	trajectory->extendKnotsTo(traj_end + 0.1 * pose_interval, Sophus::SO3d(Eigen::Quaterniond::Identity()), Eigen::Vector3d(0, 0, 0));
	printf("numKnots = %zd minTime = %.8lf maxTime = %.8lf\n", trajectory->numKnots(), trajectory->minTime(), trajectory->maxTime());
	//trajectory->print_knots();
	initTrajectory(posedata, trajectory);

	// IMU数据
	std::vector<liso::IMUData> IMUDatas;
	Eigen::Vector3d IMURotAngle = Eigen::Vector3d(180, 0, 90);
	std::string	imu_file = "C:\\Users\\Junlong\\Desktop\\2023-01-05-01\\363880-363990\\ROVE_20230105_01_STIM300.imr.txt";
	liso::readIMUData(imu_file, IMUDatas, IMURotAngle*D2R);

	// 计算残差
	FILE*	fout = fopen("imu_residual.txt", "w");
	for (size_t i = 0; i < IMUDatas.size(); i++)
	{
		if (IMUDatas[i].timestamp < traj_start)
			continue;
		if (IMUDatas[i].timestamp > traj_end)
			break;

		if (fabs(std::remainder(IMUDatas[i].timestamp, 1.0)) < 0.003)
		{
			Eigen::Vector3d acc_residual = trajectory->accelResidual(IMUDatas[i].timestamp, IMUDatas[i].accel);
			Eigen::Vector3d gyr_residual = trajectory->gyroResidual(IMUDatas[i].timestamp, IMUDatas[i].gyro);
			fprintf(fout, "%10.3lf %14.5lf %14.5lf %14.5lf %14.5lf %14.5lf %14.5lf\n", IMUDatas[i].timestamp,
				acc_residual(0), acc_residual(1), acc_residual(2), gyr_residual(0), gyr_residual(1), gyr_residual(2));
		}
	}

	if (fout) { fclose(fout); fout = NULL; }
	return;
}


static void test_gnss_sins()
{
	const double knot_interval = 0.1;
	double	pos_weight = 1.0;
	double	vel_weight = 5.0;
	double	rot_weight = 10.0;
	double	gyr_weight = 1.0 / 0.118;
	double	acc_weight = 1.0 / 0.664;
	double	GyrPNSD = 5.06797954e-06;
	double	AccPNSD = 9.00000000e-05;

	// 轨迹文件
	int		ref_type = 0; 
	std::string	pos_file = "363880-363990.fts";
	std::vector<liso::PoseData> posedata;
	liso::readPoseData(ref_type, pos_file, posedata, false);

	// GNSS
	std::vector<liso::GNSSPVA> GNSDatas;
	std::string	gns_file = "363880-363990.rtk";
	liso::readGNSSPVAData(gns_file, GNSDatas, false);

	// IMU数据
	std::vector<liso::IMUData> IMUDatas;
	Eigen::Vector3d IMURotAngle = Eigen::Vector3d(180, 0, 90);
	std::string	imu_file = "363880-363990.imu";
	liso::readIMUData(imu_file, IMUDatas, IMURotAngle*D2R);

	double traj_start = posedata[0].timestamp;
	double traj_end = posedata.back().timestamp;
	double traj_timespan = traj_end - traj_start;
	const double segment_interval = 5.0; // 分段时间间隔,主要处理惯导的零偏

	// 重新调整结束时间
	int    segment_number = (int)(traj_timespan / segment_interval);
	double	remain_time = fmod(traj_timespan, segment_interval);
	if (remain_time > 0)
	{
		// 删除掉多余的部分
		traj_end = traj_start + segment_number * segment_interval;
		traj_timespan = traj_end - traj_start;
	}

	// 分段参数 重力方向/陀螺零偏/加计零偏/激光时间偏差
	liso::CalibParamManager::Ptr calib_param = std::make_shared<liso::CalibParamManager>();
	calib_param->traj_starttime = traj_start;
	calib_param->traj_endtime = traj_end;
	calib_param->segment_interval = segment_interval;
	calib_param->segment_number = segment_number;
	for (int i = 0; i < segment_number; i++)
	{
		std::pair<double, double> segment_time;
		segment_time.first = traj_start + i * segment_interval;
		segment_time.second = traj_start + (i + 1) * segment_interval;
		calib_param->segment_timestamp.push_back(segment_time);

		liso::SegmentCalibParam	segment_param;
		segment_param.starttime = traj_start + i * segment_interval; 
		segment_param.endtime = traj_start + (i + 1) * segment_interval;
		calib_param->segment_param.emplace_back(segment_param);
		fprintf(stdout, "segment %4d %10.3lf %10.3lf\n", i, segment_time.first, segment_time.second);
	}

	
	std::shared_ptr<liso::Trajectory> trajectory = std::make_shared<liso::Trajectory>(knot_interval, traj_start, 0);
	trajectory->extendKnotsTo(traj_end + 0.1 * knot_interval, Sophus::SO3d(Eigen::Quaterniond::Identity()), Eigen::Vector3d(0, 0, 0));
	trajectory->SetCalibParam(calib_param);

	// init
	initTrajectory(posedata, trajectory);
	initIMUAccelBias(IMUDatas, trajectory, calib_param);
	initIMUGyrosBias(IMUDatas, trajectory, calib_param);

	liso::TrajectoryEstimatorOptions options;
	options.lock_trajP = false;
	options.lock_trajR = false;
	options.lock_ab = false;
	options.lock_wb = false;
	options.lock_g = true;
	options.local_traj = false;
	options.lock_tlb = true;
	options.lock_Rlb = true;


	liso::TrajectoryEstimator estimator(trajectory, calib_param, options);

	// GNSS
	Eigen::Vector3d  blarm(0.500, -0.282, -0.057);
	for (int i = 0; i < GNSDatas.size(); i++)
	{
		estimator.addGNSSPosMeasurement(GNSDatas[i], blarm, pos_weight);
		estimator.addGNSSVelMeasurement(GNSDatas[i], blarm, vel_weight);
	}

	// IMU
	for (int i = 0; i < IMUDatas.size(); i++)
	{
		estimator.AddIMUMeasurement(IMUDatas[i], gyr_weight, acc_weight);
	}
	estimator.AddGyroBiasBetweenConstraint(GyrPNSD);
	estimator.AddAccelBiasBetweenConstraint(AccPNSD);

	ceres::Solver::Summary summary = estimator.Solve(50, true);

	const std::string outPose = pos_file + ".out";
	FILE*	fout = fopen(outPose.c_str(), "w");
	for (int i = 0; i < posedata.size(); i++)
	{
		liso::SE3d pose = estimator.trajectory_->pose(posedata[i].timestamp);
		Eigen::Vector3d vel = estimator.trajectory_->velocity(posedata[i].timestamp);

		double	Azimuth[3];
		Eigen::Vector3d pos = pose.translation();
		Eigen::Matrix3d Rbe = pose.unit_quaternion().toRotationMatrix();
		liso::PoseData2Attitude(pose, Azimuth, true);
		fprintf(fout ? fout : stdout, "%d %10.3lf %14.5lf %14.5lf %14.5lf %9.6lf %9.6lf %9.6lf %9.6lf %9.6lf %9.6lf\n", int(posedata[i].timestamp / 604800.0), fmod(posedata[i].timestamp, 604800.0),
			pos(0), pos(1), pos(2), vel(0), vel(1), vel(2), Azimuth[0], Azimuth[1], Azimuth[2]);
	}
	if (fout) { fclose(fout); fout = NULL; }


	return;
}

static void test_gnss_sins_lidar()
{
	const double knot_interval = 0.1;
	double	pos_weight = 1.0;
	double	vel_weight = 5.0;
	double	rot_weight = 10.0;
	double	gyr_weight = 1.0 / 0.118;
	double	acc_weight = 1.0 / 0.664;
	double	GyrPNSD = 5.06797954e-06;
	double	AccPNSD = 9.00000000e-05;

	// 轨迹文件
	int		ref_type = 0;
	std::string	pos_file = "363880-363990.fts";
	std::vector<liso::PoseData> posedata;
	liso::readPoseData(ref_type, pos_file, posedata, false);

	// GNSS
	std::vector<liso::GNSSPVA> GNSDatas;
	std::string	gns_file = "363880-363990.rtk";
	liso::readGNSSPVAData(gns_file, GNSDatas, false);

	// IMU数据
	std::vector<liso::IMUData> IMUDatas;
	Eigen::Vector3d IMURotAngle = Eigen::Vector3d(180, 0, 90);
	std::string	imu_file = "363880-363990.imu";
	liso::readIMUData(imu_file, IMUDatas, IMURotAngle*D2R);

	double traj_start = posedata[0].timestamp;
	double traj_end = posedata.back().timestamp;
	double traj_timespan = traj_end - traj_start;
	const double segment_interval = 5.0; // 分段时间间隔,主要处理惯导的零偏

	// 重新调整结束时间
	int    segment_number = (int)(traj_timespan / segment_interval);
	double	remain_time = fmod(traj_timespan, segment_interval);
	if (remain_time > 0)
	{
		// 删除掉多余的部分
		traj_end = traj_start + segment_number * segment_interval;
		traj_timespan = traj_end - traj_start;
	}

	// 分段参数 重力方向/陀螺零偏/加计零偏/激光时间偏差
	liso::CalibParamManager::Ptr calib_param = std::make_shared<liso::CalibParamManager>();
	calib_param->traj_starttime = traj_start;
	calib_param->traj_endtime = traj_end;
	calib_param->segment_interval = segment_interval;
	calib_param->segment_number = segment_number;
	for (int i = 0; i < segment_number; i++)
	{
		std::pair<double, double> segment_time;
		segment_time.first = traj_start + i * segment_interval;
		segment_time.second = traj_start + (i + 1) * segment_interval;
		calib_param->segment_timestamp.push_back(segment_time);

		liso::SegmentCalibParam	segment_param;
		segment_param.starttime = traj_start + i * segment_interval;
		segment_param.endtime = traj_start + (i + 1) * segment_interval;
		calib_param->segment_param.emplace_back(segment_param);
		fprintf(stdout, "segment %4d %10.3lf %10.3lf\n", i, segment_time.first, segment_time.second);
	}

	
	std::shared_ptr<liso::Trajectory> trajectory = std::make_shared<liso::Trajectory>(knot_interval, traj_start, 0);
	trajectory->extendKnotsTo(traj_end + 0.1 * knot_interval, Sophus::SO3d(Eigen::Quaterniond::Identity()), Eigen::Vector3d(0, 0, 0));
	trajectory->SetCalibParam(calib_param);

	// init
	initTrajectory(posedata, trajectory);
	initIMUAccelBias(IMUDatas, trajectory, calib_param);
	initIMUGyrosBias(IMUDatas, trajectory, calib_param);

	Eigen::Matrix3d Rlb = Eigen::Matrix3d::Identity();
	Eigen::Vector3d tlb = Eigen::Vector3d(0.000, 0.000, 0.000);
	Rlb << -0.99999536609835638, -0.00140092145000101, 0.00270281351655828,
		0.00139808552502202, -0.99999847053077795, -0.00105085344816699,
		0.00270428154582451, -0.00104706981416111, 0.99999579524422266;
	tlb = Eigen::Vector3d(0.000, 0.000, 0.17445);

	lpostk::FrameContainer frameContainer;
	frameContainer.open_poles_filelist("filelst-pole.lst");
	frameContainer.open_planes_filelist("filelst-plane.lst");

	// 构建共视关系
	frameContainer.set_frame_pose(trajectory, Rlb, tlb);
	frameContainer.TrackPoint();
	frameContainer.TrackLine();
	frameContainer.TrackPole();
	frameContainer.TrackPlane();

	liso::TrajectoryEstimatorOptions options;
	options.lock_trajP = false;
	options.lock_trajR = false;
	options.lock_ab = false;
	options.lock_wb = false;
	options.lock_g = true;
	options.local_traj = false;
	options.lock_tlb = true;
	options.lock_Rlb = true;

	liso::TrajectoryEstimator estimator(trajectory, calib_param, options);

	// GNSS
	Eigen::Vector3d  blarm(0.500, -0.282, -0.057);
	for (int i = 0; i < GNSDatas.size(); i++)
	{
		estimator.addGNSSPosMeasurement(GNSDatas[i], blarm, pos_weight);
		estimator.addGNSSVelMeasurement(GNSDatas[i], blarm, vel_weight);
	}

	// IMU
	for (int i = 0; i < IMUDatas.size(); i++)
	{
		estimator.AddIMUMeasurement(IMUDatas[i], gyr_weight, acc_weight);
	}
	estimator.AddGyroBiasBetweenConstraint(GyrPNSD);
	estimator.AddAccelBiasBetweenConstraint(AccPNSD);

	// LiDAR
	estimator.pFrameContainer = &frameContainer;
	estimator.AddPointMeasurements(LiDMeasType::LOCALTOLOCAL, 0.10);
	estimator.AddLineMeasurements(LiDMeasType::LOCALTOLOCAL, 0.10, 5 * D2R);
	estimator.AddPoleMeasurements(LiDMeasType::LOCALTOLOCAL, 0.10, 5 * D2R);
	estimator.AddPlaneMeasurements(LiDMeasType::LOCALTOLOCAL, 0.10, 5 * D2R);

	ceres::Solver::Summary summary = estimator.Solve(50, true);

	const std::string outPose = pos_file + ".out";
	FILE*	fout = fopen(outPose.c_str(), "w");
	for (int i = 0; i < posedata.size(); i++)
	{
		liso::SE3d pose = estimator.trajectory_->pose(posedata[i].timestamp);
		Eigen::Vector3d vel = estimator.trajectory_->velocity(posedata[i].timestamp);

		double	Azimuth[3];
		Eigen::Vector3d pos = pose.translation();
		Eigen::Matrix3d Rbe = pose.unit_quaternion().toRotationMatrix();
		liso::PoseData2Attitude(pose, Azimuth, true);
		fprintf(fout ? fout : stdout, "%d %10.3lf %14.5lf %14.5lf %14.5lf %9.6lf %9.6lf %9.6lf %9.6lf %9.6lf %9.6lf\n", int(posedata[i].timestamp / 604800.0), fmod(posedata[i].timestamp, 604800.0),
			pos(0), pos(1), pos(2), vel(0), vel(1), vel(2), Azimuth[0], Azimuth[1], Azimuth[2]);
	}
	if (fout) { fclose(fout); fout = NULL; }


	return;
}


int main(int argc, char *argv[])
{	
  test_gnss_sins_lidar();
	return 0;
}