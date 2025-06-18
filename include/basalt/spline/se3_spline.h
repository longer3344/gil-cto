/**
BSD 3-Clause License

This file is part of the Basalt project.
https://gitlab.com/VladyslavUsenko/basalt-headers.git

Copyright (c) 2019, Vladyslav Usenko and Nikolaus Demmel.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

@file
@brief Uniform B-spline for SE(3)
*/

#pragma once

#include <basalt/spline/rd_spline.h>
#include <basalt/spline/so3_spline.h>
#include <basalt/utils/assert.h>

#include <basalt/spline/calib_bias.hpp>
#include <basalt/spline/spline_segment.h>

#include <array>
#include "CommonFunc.h"

namespace basalt {

	/// @brief Uniform B-spline for SE(3) of order N. Internally uses an SO(3) (\ref
	/// So3Spline) spline for rotation and 3D Euclidean spline (\ref RdSpline) for
	/// translation (split representaion).
	///
	/// See [[arXiv:1911.08860]](https://arxiv.org/abs/1911.08860) for more details.
	/// _N       样条阶数, 对于Bezier曲线，阶数和次数是一样的，但是B样条，阶数是次数加1
	/// _Scalar  数据类型, 一般为double
	template <int _N, typename _Scalar = double>
	class Se3Spline {
	public:
		static constexpr int N = _N;        ///< 样条阶数, 对于Bezier曲线，阶数和次数是一样的，但是B样条，阶数是次数加1
		static constexpr int DEG = _N - 1;  ///< 样条次数, 对于Bezier曲线，阶数和次数是一样的，但是B样条，阶数是次数加1

		using MatN = Eigen::Matrix<_Scalar, _N, _N>;
		using VecN = Eigen::Matrix<_Scalar, _N, 1>;
		using VecNp1 = Eigen::Matrix<_Scalar, _N + 1, 1>;

		using Vec3 = Eigen::Matrix<_Scalar, 3, 1>;
		using Vec6 = Eigen::Matrix<_Scalar, 6, 1>;
		using Vec9 = Eigen::Matrix<_Scalar, 9, 1>;
		using Vec12 = Eigen::Matrix<_Scalar, 12, 1>;

		using Mat3 = Eigen::Matrix<_Scalar, 3, 3>;
		using Mat6 = Eigen::Matrix<_Scalar, 6, 6>;

		using Mat36 = Eigen::Matrix<_Scalar, 3, 6>;
		using Mat39 = Eigen::Matrix<_Scalar, 3, 9>;
		using Mat312 = Eigen::Matrix<_Scalar, 3, 12>;

		using Matrix3Array = std::array<Mat3, N>;
		using Matrix36Array = std::array<Mat36, N>;
		using Matrix6Array = std::array<Mat6, N>;

		using SO3 = Sophus::SO3<_Scalar>;
		using SE3 = Sophus::SE3<_Scalar>;

		using PosJacobianStruct = typename RdSpline<3, N, _Scalar>::JacobianStruct;
		using SO3JacobianStruct = typename So3Spline<N, _Scalar>::JacobianStruct;

		/// @brief Struct to store the accelerometer residual Jacobian with
		/// respect to knots
		struct AccelPosSO3JacobianStruct {
			size_t start_idx;
			std::array<Mat36, N> d_val_d_knot;
		};

		/// @brief Struct to store the pose Jacobian with respect to knots
		struct PosePosSO3JacobianStruct {
			size_t start_idx;
			std::array<Mat6, N> d_val_d_knot;
		};

		/// @brief Constructor with knot interval and start time
		///
		/// @param[in] time_interval knot time interval in seconds
		/// @param[in] start_time start time of the spline in seconds
		Se3Spline(double time_interval, double start_time = 0)
			: pos_spline(time_interval, start_time),
			so3_spline(time_interval, start_time),
			dt_(time_interval) {}

		/// @brief Gererate random trajectory
		///
		/// @param[in] n number of knots to generate
		/// @param[in] static_init if true the first N knots will be the same
		/// resulting in static initial condition
		void genRandomTrajectory(int n, bool static_init = false) {
			so3_spline.genRandomTrajectory(n, static_init);
			pos_spline.genRandomTrajectory(n, static_init);
		}

		/// @brief Set the knot to particular SE(3) pose
		///
		/// @param[in] pose SE(3) pose
		/// @param[in] i index of the knot
		void setKnot(const SE3 &pose, int i) {
			so3_spline.getKnot(i) = pose.so3();
			pos_spline.getKnot(i) = pose.translation();
		}

		/// @brief Set the knot to particular Vec3 pose
		///
		/// @param[in] pos Vec3 pose
		/// @param[in] i index of the knot
		void setKnotPos(const Vec3 pos, int i) {
			pos_spline.getKnot(i) = pos;
		}

		/// @brief Set the knot to particular SO3 pose
		///
		/// @param[in] ori SO3 pose
		/// @param[in] i index of the knot
		void setKnotSO3(const SO3 ori, int i) {
			so3_spline.getKnot(i) = ori;
		}

		/// @brief Reset spline to have num_knots initialized at pose
		///
		/// @param[in] pose SE(3) pose
		/// @param[in] num_knots number of knots to initialize
		void setKnots(const SE3 &pose, int num_knots) {
			so3_spline.resize(num_knots);
			pos_spline.resize(num_knots);

			for (int i = 0; i < num_knots; i++) {
				so3_spline.getKnot(i) = pose.so3();
				pos_spline.getKnot(i) = pose.translation();
			}
		}

		/// @brief Reset spline to the knots from other spline
		///
		/// @param[in] other spline to copy knots from
		void setKnots(const Se3Spline<N, _Scalar> &other) {
			BASALT_ASSERT(other.dt_ == dt_);
			BASALT_ASSERT(other.pos_spline.getKnots().size() ==
				other.pos_spline.getKnots().size());

			size_t num_knots = other.pos_spline.getKnots().size();

			so3_spline.resize(num_knots);
			pos_spline.resize(num_knots);

			for (size_t i = 0; i < num_knots; i++) {
				so3_spline.getKnot(i) = other.so3_spline.getKnot(i);
				pos_spline.getKnot(i) = other.pos_spline.getKnot(i);
			}
		}

		/// @brief extend trajectory to time t
		///
		/// @param[in] t time
		/// @param[in] initial_so3 initial knot of so3_spline
		/// @param[in] initial_pos initial knot of pos_spline
		void extendKnotsTo(double t, const SO3 initial_so3, const Vec3 initial_pos) {
			// numKnots() < N : 开头插入N个空的数据体, 最终节点的个数为N+P,N为阶数(次数+1),P为控制点个数
			// maxTime() < t  : 为了保证最后一个时刻观测值也能使用, 需要向后扩展一个knot
			while ((numKnots() < N) || (maxTime() < t)) {
				so3_spline.knots_push_back(initial_so3);
				pos_spline.knots_push_back(initial_pos);
				//printf("numKnots = %zd minTime = %.8lf maxTime = %.8lf\n", numKnots(), pos_spline.minTime(), pos_spline.maxTime());
			}
		}


		/// @brief extend trajectory to time t
		///
		/// @param[in] t timestamp
		/// @param[in] initial_knot initial knot
		void extendKnotsTo(double timestamp, const SE3 & initial_knot) {
			// numKnots() < N : 开头插入N个空的数据体
			while ((numKnots() < N) || (maxTime() < timestamp)) {
				knots_push_back(initial_knot);
			}
		}

		/// @brief Add knot to the end of the spline
		///
		/// @param[in] knot knot to add
		inline void knots_push_back(const SE3 &knot) {
			so3_spline.knots_push_back(knot.so3());
			pos_spline.knots_push_back(knot.translation());
		}

		/// @brief Remove knot from the back of the spline
		inline void knots_pop_back() {
			so3_spline.knots_pop_back();
			pos_spline.knots_pop_back();
		}

		/// @brief Return the first knot of the spline
		///
		/// @return first knot of the spline
		inline SE3 knots_front() const {
			SE3 res(so3_spline.knots_front(), pos_spline.knots_front());
			return res;
		}

		/// @brief Remove first knot of the spline and increase the start time
		inline void knots_pop_front() {
			so3_spline.knots_pop_front();
			pos_spline.knots_pop_front();

			BASALT_ASSERT(so3_spline.minTime() == pos_spline.minTime());
			BASALT_ASSERT(so3_spline.getKnots().size() == pos_spline.getKnots().size());
		}

		/// @brief Return the last knot of the spline
		///
		/// @return last knot of the spline
		SE3 getLastKnot() {
			BASALT_ASSERT(so3_spline.getKnots().size() == pos_spline.getKnots().size());
			SE3 res(so3_spline.getKnots().back(), pos_spline.getKnots().back());
			return res;
		}

		/// @brief Return knot with index i
		///
		/// @param i index of the knot
		/// @return knot
		SE3 getKnot(size_t i) const {
			SE3 res(getKnotSO3(i), getKnotPos(i));
			return res;
		}

		/// @brief Return reference to the SO(3) knot with index i
		///
		/// @param i index of the knot
		/// @return reference to the SO(3) knot
		inline SO3 &getKnotSO3(size_t i) {
			return so3_spline.getKnot(i);
		}

		/// @brief Return const reference to the SO(3) knot with index i
		///
		/// @param i index of the knot
		/// @return const reference to the SO(3) knot
		inline const SO3 &getKnotSO3(size_t i) const {
			return so3_spline.getKnot(i);
		}

		/// @brief Return reference to the position knot with index i
		///
		/// @param i index of the knot
		/// @return reference to the position knot
		inline Vec3 &getKnotPos(size_t i) {
			return pos_spline.getKnot(i);
		}

		/// @brief Return const reference to the position knot with index i
		///
		/// @param i index of the knot
		/// @return const reference to the position knot
		inline const Vec3 &getKnotPos(size_t i) const {
			return pos_spline.getKnot(i);
		}

		/// @brief Set start time for spline
		///
		/// @param[in] start_time start time of the spline in seconds
		inline void setStartTime(double timestamp) {
			so3_spline.setStartTime(timestamp);
			pos_spline.setStartTime(timestamp);
		}

		/// @brief Apply increment to the knot
		///
		/// The incremernt vector consists of translational and rotational parts \f$
		/// [\upsilon, \omega]^T \f$. Given the current pose of the knot \f$ R \in
		/// SO(3), p \in \mathbb{R}^3\f$ the updated pose is: \f{align}{ R' &=
		/// \exp(\omega) R
		/// \\ p' &= p + \upsilon
		/// \f}
		///  The increment is consistent with \ref
		/// PoseState::applyInc.
		///
		/// @param[in] i index of the knot
		/// @param[in] inc 6x1 increment vector
		template <typename Derived>
		void applyInc(int i, const Eigen::MatrixBase<Derived> &inc) {
			EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Derived, 6);

			pos_spline.getKnot(i) += inc.template head<3>();
			so3_spline.getKnot(i) =
				SO3::exp(inc.template tail<3>()) * so3_spline.getKnot(i);
		}

		double time_interval()
		{
			return dt_;
		}

		/// @brief Maximum time represented by spline
		///
		/// @return maximum time represented by spline in nanoseconds
		/// 实际时间要小于maxTime
		double maxTime() const {
			BASALT_ASSERT_STREAM(so3_spline.maxTime() == pos_spline.maxTime(),
				"so3_spline.maxTime() " << so3_spline.maxTime()
				<< " pos_spline.maxTime() "
				<< pos_spline.maxTime());
			return pos_spline.maxTime();
		}

		/// @brief Minimum time represented by spline
		///
		/// @return minimum time represented by spline in seconds
		/// 实际时间要大于等于minTime
		double minTime() const {
			BASALT_ASSERT_STREAM(so3_spline.minTime() == pos_spline.minTime(),
				"so3_spline.minTime() " << so3_spline.minTime()
				<< " pos_spline.minTime() "
				<< pos_spline.minTime());
			return pos_spline.minTime();
		}

		double getTime(const int i) const
		{
			BASALT_ASSERT_STREAM(so3_spline.minTime() == pos_spline.minTime(),
				"so3_spline.minTime() " << so3_spline.minTime()
				<< " pos_spline.minTime() "
				<< pos_spline.minTime());

			return pos_spline.getTime(i);
		}

		/// @brief Number of knots in the spline
		size_t numKnots() const { return pos_spline.getKnots().size(); }

		/// @brief Linear acceleration in the world frame.
		///
		/// @param[in] time time to evaluate linear acceleration in seconds
		inline Vec3 transAccelWorld(double time) const {
			return pos_spline.acceleration(time);
		}

		/// @brief Linear velocity in the world frame.
		///
		/// @param[in] time time to evaluate linear velocity in seconds
		inline Vec3 transVelWorld(double time) const {
			return pos_spline.velocity(time);
		}

		/// @brief Rotational velocity in the body frame.
		///
		/// @param[in] time time to evaluate rotational velocity in seconds
		inline Vec3 rotVelBody(double time) const {
			return so3_spline.velocityBody(time);
		}

		/// @brief Evaluate pose.
		///
		/// @param[in] time time to evaluate pose in seconds
		/// @return SE(3) pose at time
		SE3 pose(double time) const {
			SE3 res;
			res.so3() = so3_spline.evaluate(time);
			res.translation() = pos_spline.evaluate(time);

			return res;
		}

		/// @brief Evaluate pose and compute Jacobian.
		///
		/// @param[in] time time to evaluate pose inseconds
		/// @param[out] J Jacobian of the pose with respect to knots
		/// @return SE(3) pose at time
		SE3 pose(double time, PosePosSO3JacobianStruct *J) const {
			SE3 res;

			typename So3Spline<_N, _Scalar>::JacobianStruct Jr;
			typename RdSpline<3, N, _Scalar>::JacobianStruct Jp;

			res.so3() = so3_spline.evaluate(time, &Jr);
			res.translation() = pos_spline.evaluate(time, &Jp);

			if (J) {
				Eigen::Matrix3d RT = res.so3().inverse().matrix();

				J->start_idx = Jr.start_idx;
				for (int i = 0; i < N; i++) {
					J->d_val_d_knot[i].setZero();
					J->d_val_d_knot[i].template topLeftCorner<3, 3>() = RT * Jp.d_val_d_knot[i];     // 姿态
					J->d_val_d_knot[i].template bottomRightCorner<3, 3>() = RT * Jr.d_val_d_knot[i]; // 位置
				}
			}

			return res;
		}

		/// @brief Evaluate pose and compute time Jacobian.
		///
		/// @param[in] time time to evaluate pose in seconds
		/// @param[out] J Jacobian of the pose with time
		void d_pose_d_t(double time, Vec6 &J) const {
			J.template head<3>() = so3_spline.evaluate(time).inverse() * transVelWorld(time);
			J.template tail<3>() = rotVelBody(time);
		}

		/// @brief 计算全局系下的位置
		///
		/// @param[in] time         当前时间
		/// @return Vec3
		Vec3 position(double time) const
		{
			Vec3 pos_world = pos_spline.evaluate<0>(time);
			return pos_world;
		}

		/// @brief 计算全局系下的速度
		///
		/// @param[in] time         当前时间
		/// @return Vec3
		Vec3 velocity(double time) const
		{
			Vec3 vel_world = pos_spline.velocity(time);
			return vel_world;
		}

		/// @brief 计算惯导系下的速度
		///
		/// @param[in] time         当前时间
		/// @return Vec3
		Vec3 velocity_body(double time) const {
			Sophus::SO3d R = so3_spline.evaluate(time);
			Vec3 vel_world = pos_spline.velocity(time);
			return R.inverse() * vel_world;
		}

		/// @brief 计算全局系下的姿态
		///
		/// @param[in] time         当前时间
		/// @return SO3d
		Sophus::SO3d orientation(double time)
		{
			return so3_spline.evaluate(time);
		}

		/// @brief 计算全局系下的加速度
		///
		/// @param[in] time         当前时间
		/// @return Vec3
		Vec3 accel(double time) const
		{
			Eigen::Vector3d accel_world = pos_spline.acceleration(time);
			return accel_world;
		}

		/// @brief 计算惯导系下的加速度
		///
		/// @param[in] time         当前时间
		/// @param[in] gravity      重力加速度
		/// @return Vec3
		Vec3 accel_body(double time) const
		{
			// 导出位置/速度/加速度/姿态
			Eigen::Vector3d pos_world = pos_spline.evaluate<0>(time);
			Eigen::Vector3d vel_world = pos_spline.evaluate<1>(time);
			Eigen::Vector3d acc_world = pos_spline.evaluate<2>(time);
			Sophus::SO3d R_IkToG = so3_spline.evaluate(time); // 姿态

			// 重力
			Eigen::Vector3d LLH = XYZ2LLH(pos_world);
			Eigen::Vector3d gravity = getGravityECEF(LLH);

			if (time == 368970.0)
			{
				std::cout << time << std::endl;
				std::cout << "gravity " << std::endl << gravity << std::endl;
				std::cout << "accel_world " << std::endl << acc_world << std::endl;
				std::cout << "accel_world - gravity " << std::endl << acc_world - gravity << std::endl;

				Eigen::Quaterniond qqq = R_IkToG.unit_quaternion();
				std::cout << qqq.x() << qqq.y() << qqq.z() << qqq.w() << std::endl;
			}

			// 科氏力
			Eigen::Vector3d Wie(0.0, 0.0, 7.2921151467E-5);
			Eigen::Vector3d coriolis = 2 * Wie.cross(vel_world);

			Eigen::Vector3d acc_body = R_IkToG.inverse() * (acc_world + coriolis - gravity); // 注意g的方向
			return acc_body;
		}

		/// @brief 计算惯导系下的加速度
		///
		/// @param[in] time         当前时间
		/// @param[in] gravity      重力加速度
		/// @return Vec3
		Vec3 accel_body(double time, Sophus::SO3d &R_IkToG) const
		{
			// 导出位置/速度/加速度/姿态
			Eigen::Vector3d pos_world = pos_spline.evaluate<0>(time);
			Eigen::Vector3d vel_world = pos_spline.evaluate<1>(time);
			Eigen::Vector3d acc_world = pos_spline.evaluate<2>(time);

			// 重力
			Eigen::Vector3d LLH = XYZ2LLH(pos_world);
			Eigen::Vector3d gravity = getGravityECEF(LLH);

			if (time == 368970.0)
			{
				std::cout << time << std::endl;
				std::cout << "gravity " << std::endl << gravity << std::endl;
				std::cout << "accel_world " << std::endl << acc_world << std::endl;
				std::cout << "accel_world - gravity " << std::endl << acc_world - gravity << std::endl;

				Eigen::Quaterniond qqq = R_IkToG.unit_quaternion();
				std::cout << qqq.x() << qqq.y() << qqq.z() << qqq.w() << std::endl;
			}

			// 科氏力
			Eigen::Vector3d Wie(0.0, 0.0, 7.2921151467E-5);
			Eigen::Vector3d coriolis = 2 * Wie.cross(vel_world);

			Eigen::Vector3d acc_body = R_IkToG.inverse() * (acc_world + coriolis - gravity); // 注意g的方向
			return acc_body;
		}

		/// @brief 计算惯导系下的加速度
		///
		/// @param[in] time         当前时间
		/// @param[in] gravity      重力加速度
		/// @return Vec3
		Vec3 accel_body(double time, const Eigen::Vector3d &gravity) const
		{
			// 导出位置/速度/加速度/姿态
			Eigen::Vector3d pos_world = pos_spline.evaluate<0>(time);
			Eigen::Vector3d vel_world = pos_spline.evaluate<1>(time);
			Eigen::Vector3d acc_world = pos_spline.evaluate<2>(time);
			Sophus::SO3d R_IkToG = so3_spline.evaluate(time); // 姿态

			if (time == 368970.0)
			{
				std::cout << time << std::endl;
				std::cout << "gravity " << std::endl << gravity << std::endl;
				std::cout << "accel_world " << std::endl << acc_world << std::endl;
				std::cout << "accel_world - gravity " << std::endl << acc_world - gravity << std::endl;

				Eigen::Quaterniond qqq = R_IkToG.unit_quaternion();
				std::cout << qqq.x() << qqq.y() << qqq.z() << qqq.w() << std::endl;
			}

			// 科氏力
			Eigen::Vector3d Wie(0.0, 0.0, 7.2921151467E-5);
			Eigen::Vector3d coriolis = 2 * Wie.cross(vel_world);

			Eigen::Vector3d acc_body = R_IkToG.inverse() * (acc_world + coriolis - gravity); // 注意g的方向
			return acc_body;
		}

		/// @brief 计算惯导系下的加速度
		///
		/// @param[in] time         当前时间
		/// @param[in] R_IkToG      姿态
		/// @param[in] gravity      重力加速度
		/// @return Vec3
		Vec3 accel_body(double time, Sophus::SO3d &R_IkToG, const Eigen::Vector3d &gravity) const
		{
			// 导出位置/速度/加速度
			Eigen::Vector3d pos_world = pos_spline.evaluate<0>(time);
			Eigen::Vector3d vel_world = pos_spline.evaluate<1>(time);
			Eigen::Vector3d acc_world = pos_spline.evaluate<2>(time);

			if (time == 368970.0)
			{
				std::cout << time << std::endl;
				std::cout << "gravity " << std::endl << gravity << std::endl;
				std::cout << "accel_world " << std::endl << acc_world << std::endl;
				std::cout << "accel_world - gravity " << std::endl << acc_world - gravity << std::endl;

				Eigen::Quaterniond qqq = R_IkToG.unit_quaternion();
				std::cout << qqq.x() << qqq.y() << qqq.z() << qqq.w() << std::endl;
			}

			// 科氏力
			Eigen::Vector3d Wie(0.0, 0.0, 7.2921151467E-5);
			Eigen::Vector3d coriolis = 2 * Wie.cross(vel_world);

			Eigen::Vector3d acc_body = R_IkToG.inverse() * (acc_world + coriolis - gravity); // 注意g的方向
			return acc_body;
		}

		Vec3 wbb(double time) const {
			return so3_spline.velocityBody(time);
		}

		/// @brief 计算陀螺残差
		///
		/// @param[in] time          当前时间
		/// @param[in] measurement   陀螺原始观测
		/// @return Vec3
		Vec3 gyroResidual(double time, const Vec3 &measurement, const Eigen::Vector3d &bias = Eigen::Vector3d(0.0, 0.0, 0.0)) const
		{
			Eigen::Vector3d Wie(0.0, 0.0, 7.2921151467E-5);     // 地球自传
			Sophus::SO3d rot_world = so3_spline.evaluate(time); // 姿态
			return so3_spline.velocityBody(time) + rot_world * Wie - measurement + bias;
		}

		/// @brief 计算陀螺残差
		///
		/// @param[in] time           当前时间
		/// @param[in] measurement    陀螺原始观测
		/// @param[in] gyro_bias_full 陀螺零偏/比例因子/非正交
		/// @return Vec3
		Vec3 gyroResidual(double time, const Vec3 &measurement, const CalibGyroBias<_Scalar> &gyro_bias_full) const
		{
			return so3_spline.velocityBody(time) - gyro_bias_full.getCalibrated(measurement);
		}

		/// @brief Evaluate gyroscope residual and compute Jacobians.
		///
		/// @param[in] time time of the measurement
		/// @param[in] measurement gyroscope measurement
		/// @param[in] gyro_bias_full gyroscope calibration
		/// @param[out] J_knots Jacobian with respect to SO(3) spline knots
		/// @param[out] J_bias Jacobian with respect to gyroscope calibration
		/// @return gyroscope residual
		Vec3 gyroResidual(double time, const Vec3 &measurement, const CalibGyroBias<_Scalar> &gyro_bias_full,
			SO3JacobianStruct *J_knots, Mat312 *J_bias = nullptr) const
		{
			if (J_bias)
			{
				J_bias->setZero();
				J_bias->template block<3, 3>(0, 0).diagonal().array() = 1.0; // bx/by/bz
				J_bias->template block<3, 3>(0, 3).diagonal().array() = -measurement[0]; // s1/s2/s3
				J_bias->template block<3, 3>(0, 6).diagonal().array() = -measurement[1]; // s4/s5/s6
				J_bias->template block<3, 3>(0, 9).diagonal().array() = -measurement[2]; // s7/s8/s9
			}

			///  | s1+1  s4    s7  |   |wx|   |bx|
			///  |  s2  s5+1   s8  | * |wy| - |by|
			///  |  s3   s6   s9+1 |   |wz|   |bz|
			return so3_spline.velocityBody(time, J_knots) - gyro_bias_full.getCalibrated(measurement);
		}

		/// @brief 计算加计残差
		///
		/// @param[in] time            当前时间
		/// @param[in] measurement     加计原始数据
		/// @return accelerometer residual
		Vec3 accelResidual(double time, const Eigen::Vector3d &measurement, const Eigen::Vector3d &bias = Eigen::Vector3d(0.0, 0.0, 0.0)) const
		{
			Eigen::Vector3d acc = accel_body(time);
			return acc - measurement + bias;
		}

		/// @brief 计算加计残差
		///
		/// @param[in] time            当前时间
		/// @param[in] measurement     加计原始数据
		/// @param[in] g               重力加速度
		/// @return accelerometer residual
		Vec3 accelResidualWithGravity(double time, const Eigen::Vector3d &measurement, const Eigen::Vector3d &g, 
			const Eigen::Vector3d &bias = Eigen::Vector3d(0.0, 0.0, 0.0)) const
		{
			Eigen::Vector3d acc = accel_body(time, g);
			return acc - measurement + bias;
		}

		/// @brief 计算加计残差
		///
		/// @param[in] time            当前时间
		/// @param[in] measurement     加计原始数据
		/// @param[in] accel_bias_full 加计零偏/比例因子/非正交
		/// @param[in] g               重力加速度
		/// @return accelerometer residual
		Vec3 accelResidual(double time, const Eigen::Vector3d &measurement,
			const CalibAccelBias<_Scalar> &accel_bias_full, const Eigen::Vector3d &g) const
		{
			// 导出位置/速度/加速度/姿态
			Eigen::Vector3d vel_world = pos_spline.velocity(time);
			Eigen::Vector3d accel_world = pos_spline.acceleration(time);
			Sophus::SO3d R = so3_spline.evaluate(time);

			// 科氏力
			Eigen::Vector3d Wie(0.0, 0.0, 7.2921151467E-5);
			Eigen::Vector3d coriolis = 2 * Wie.cross(vel_world);

			return R.inverse() * (accel_world + coriolis - g) - accel_bias_full.getCalibrated(measurement);
		}

		/// @brief Evaluate accelerometer residual and Jacobians.
		///
		/// @param[in] time time of the measurement
		/// @param[in] measurement accelerometer measurement
		/// @param[in] accel_bias_full accelerometer calibration
		/// @param[in] g gravity
		/// @param[out] J_knots Jacobian with respect to spline knots
		/// @param[out] J_bias Jacobian with respect to accelerometer calibration
		/// @param[out] J_g Jacobian with respect to gravity
		/// @return accelerometer residual
		Vec3 accelResidual(double time, const Vec3 &measurement,
			const CalibAccelBias<_Scalar> &accel_bias_full,
			const Vec3 &g, AccelPosSO3JacobianStruct *J_knots,
			Mat39 *J_bias = nullptr, Mat3 *J_g = nullptr) const
		{
			typename So3Spline<_N, _Scalar>::JacobianStruct Jr;
			typename RdSpline<3, N, _Scalar>::JacobianStruct Jv;
			typename RdSpline<3, N, _Scalar>::JacobianStruct Jp;

			Sophus::SO3d R = so3_spline.evaluate(time, &Jr);
			Eigen::Vector3d vel_world = pos_spline.velocity(time, Jv);
			Eigen::Vector3d accel_world = pos_spline.acceleration(time, &Jp);

			Eigen::Matrix3d RT = R.inverse().matrix();
			Eigen::Matrix3d tmp = RT * Sophus::SO3d::hat(accel_world + g);

			BASALT_ASSERT_STREAM(
				Jr.start_idx == Jp.start_idx,
				"Jr.start_idx " << Jr.start_idx << " Jp.start_idx " << Jp.start_idx);

			BASALT_ASSERT_STREAM(
				so3_spline.getKnots().size() == pos_spline.getKnots().size(),
				"so3_spline.getKnots().size() " << so3_spline.getKnots().size()
				<< " pos_spline.getKnots().size() " << pos_spline.getKnots().size());

			J_knots->start_idx = Jp.start_idx;
			for (int i = 0; i < N; i++) {
				J_knots->d_val_d_knot[i].template topLeftCorner<3, 3>() = RT * Jp.d_val_d_knot[i];
				J_knots->d_val_d_knot[i].template bottomRightCorner<3, 3>() = tmp * Jr.d_val_d_knot[i];
			}

			// 加计零偏/比例因子/非正交
			if (J_bias)
			{
				J_bias->setZero();
				J_bias->template block<3, 3>(0, 0).diagonal().array() = 1.0; // bx/by/bz
				J_bias->template block<3, 3>(0, 3).diagonal().array() = -measurement[0]; // s1/s2/s3
				(*J_bias)(1, 6) = -measurement[1]; // s4
				(*J_bias)(2, 7) = -measurement[1]; // s5
				(*J_bias)(2, 8) = -measurement[2]; // s6
			}

			// 重力加速度
			if (J_g) (*J_g) = RT;

			// 科氏力
			Eigen::Vector3d Wie(0.0, 0.0, 7.2921151467E-5);
			Vec3 coriolis = 2 * Wie.cross(vel_world);

			///  | s1+1   0    0  |   |ax|   |bx|
			///  |  s2  s4+1   0  | * |ay| - |by|
			///  |  s3   s5  s6+1 |   |az|   |bz|
			Vec3 res = RT * (accel_world - g) - accel_bias_full.getCalibrated(measurement);
			return res;
		}

		/// @brief Evaluate position residual.
		///
		/// @param[in] time time of the measurement
		/// @param[in] measured_position position measurement
		/// @param[out] Jp if not nullptr, Jacobian with respect to knos of the
		/// position spline
		/// @return position residual
		Sophus::Vector3d positionResidual(double time, const Vec3 &measured_position, PosJacobianStruct *Jp = nullptr) const
		{
			return pos_spline.evaluate(time, Jp) - measured_position;
		}

		/// @brief Evaluate orientation residual.
		///
		/// @param[in] time time of the measurement
		/// @param[in] measured_orientation orientation measurement
		/// @param[out] Jr if not nullptr, Jacobian with respect to knos of the
		/// SO(3) spline
		/// @return orientation residual
		Sophus::Vector3d orientationResidual(double time, const SO3 &measured_orientation, SO3JacobianStruct *Jr = nullptr) const
		{
			Sophus::Vector3d res = (so3_spline.evaluate(time, Jr) * measured_orientation.inverse()).log();

			if (Jr) {
				Eigen::Matrix3d Jrot;
				Sophus::leftJacobianSO3(res, Jrot);

				for (int i = 0; i < N; i++) {
					Jr->d_val_d_knot[i] = Jrot * Jr->d_val_d_knot[i];
				}
			}

			return res;
		}

		/// @brief Print knots for debugging.
		inline void print_knots() const
		{
			for (size_t i = 0; i < pos_spline.getKnots().size(); i++) {
				//std::cout << i << ": p:" << pos_spline.getKnot(i).transpose() << " q: "
				//          << so3_spline.getKnot(i).unit_quaternion().coeffs().transpose()
				//          << std::endl;

				double	time = getTime(i);
				Eigen::Vector3d pos = pos_spline.getKnot(i);
				Eigen::Quaterniond quat = so3_spline.getKnot(i).unit_quaternion();
				fprintf(stdout, "%10zd %15.8lf %14.5lf %14.5lf %14.5lf %9.6lf %9.6lf %9.6lf %9.6lf\n", i, time, pos(0), pos(1), pos(2), quat.x(), quat.y(), quat.z(), quat.w());
			}
		}

		/// @brief Print position knots for debugging.
		inline void print_pos_knots() const {
			for (size_t i = 0; i < pos_spline.getKnots().size(); i++) {
				std::cout << pos_spline.getKnot(i).transpose() << std::endl;
			}
		}

		/// @brief Knot time interval in nanoseconds.
		inline double getDt() const { return dt_; }


		/// @brief 计算当前时间的样条节点位置
		///
		/// @param[in] timestamp	当前时间
		/// @return (归一化时间,开始节点标号)
		std::pair<double, size_t> computeTIndex(double timestamp) const {
			return pos_spline.computeTIndex(timestamp);
		}

		/// @brief 获取时间范围内所有的节点
		///
		/// @param[in]  timestamp    开始时间和结束时间
		/// @param[out] spline_meta  时间范围内所有的节点
		/// @return (归一化时间,开始节点标号)
		void CaculateSplineMeta(time_init_t times, SplineMeta<_N>& spline_meta) const {
			double master_dt = getDt();   // 样条节点的时间间隔
			double master_t0 = minTime(); // 样条节点的开始时间
			size_t current_segment_start = 0;
			size_t current_segment_end = 0; // Negative signals no segment created yet

			// Times are guaranteed to be sorted correctly and t2 >= t1
			for (auto tt : times) {
				std::pair<double, size_t> ui_1, ui_2;
				ui_1 = pos_spline.computeTIndex(tt.first);  // 开始时间的样条节点位置和归一化时间
				ui_2 = pos_spline.computeTIndex(tt.second); // 结束时间的样条节点位置和归一化时间

				// 开始和结束时间的样条节点位置
				size_t i1 = ui_1.second;
				size_t i2 = ui_2.second;

				// Create new segment, or extend the current one
				if (spline_meta.segments.empty() || i1 > current_segment_end)
				{
					double segment_t0 = master_t0 + master_dt * double(i1);
					spline_meta.segments.push_back(SplineSegmentMeta<_N>(segment_t0, master_dt));
					current_segment_start = i1;
				}
				else {
					i1 = current_segment_end + 1;
				}

				// 节点个数, 需要扩展N个, N为样条阶数, 对于Bezier曲线，阶数和次数是一样的，但是B样条，阶数是次数加1
				auto& current_segment_meta = spline_meta.segments.back(); // 最后一个时间段
				for (size_t i = i1; i < (i2 + N); ++i) {
					current_segment_meta.n += 1;
				}

				current_segment_end = current_segment_start + current_segment_meta.n - 1;
			} // for times

		}

		void ShowMinMaxTime()
		{
			printf("minTime = %.3lf maxTime = %.3lf\n", pos_spline.minTime(), pos_spline.maxTime());
		}

		EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	public:
		double dt_;                           ///< 节点的时间间隔
		RdSpline<3, _N, _Scalar> pos_spline;  ///< 位置样条
		So3Spline<_N, _Scalar> so3_spline;    ///< 旋转样条

	};

}  // namespace basalt
