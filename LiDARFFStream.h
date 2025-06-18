#ifndef _LIDARFFSTREAM_H_H_
#define _LIDARFFSTREAM_H_H_


#include <string>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <algorithm>

#include "CommonFunc.h"
#include "trajectory/se3_trajectory.h"

namespace lpostk
{
	class Frame;
	class MapPoint;

	/** @brief  Subtracts two matrices : M3 = M1 - M2 */
	inline void M31_M31(const double M1[3], const double M2[3], double M3[3])
	{
		for (int i = 0; i < 3; i++) M3[i] = M1[i] - M2[i];
	}

	/** @brief  Assigns matrix : M2 = M */
	inline void M31EQU(const double M[3], double M2[3])
	{
		for (int i = 0; i < 3; i++) M2[i] = M[i];
	}


	/** @brief  Assigns matrix : M2 = -M */
	inline void M31EQU_1(const double M[3], double M2[3])
	{
		for (int i = 0; i < 3; i++) M2[i] = -M[i];
	}

	/** @brief  Scale matrix : M2 = a*M */
	inline void M31Scale(double a, double M[3])
	{
		for (int i = 0; i < 3; i++) M[i] *= a;
	}

	/** @brief  Cross two vectors : v3 = v1 X v2 */
	inline void CrossM3(const double v1[3], const double v2[3], double v3[3])
	{
		v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
		v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
		v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
	}

	/** @brief Assigns matrix : M2 = M */
	inline void M33EQU(const double M[9], double M2[9])
	{
		for (int i = 0; i < 9; i++) M2[i] = M[i];
	}

	/** @brief  Multiplys two matrices : M31_ = M33 * M31 */
	inline void M33XM31(const double M33[9], const double M31[3], double M31_[3])
	{
		M31_[0] = M33[0] * M31[0] + M33[1] * M31[1] + M33[2] * M31[2];
		M31_[1] = M33[3] * M31[0] + M33[4] * M31[1] + M33[5] * M31[2];
		M31_[2] = M33[6] * M31[0] + M33[7] * M31[1] + M33[8] * M31[2];
	}


	/** @brief  Multiplys two matrices : M31_ = M33 * M31 */
	inline void M33XM31(const double M33[9], const double M31[3], float M31_[3])
	{
		M31_[0] = float(M33[0] * M31[0] + M33[1] * M31[1] + M33[2] * M31[2]);
		M31_[1] = float(M33[3] * M31[0] + M33[4] * M31[1] + M33[5] * M31[2]);
		M31_[2] = float(M33[6] * M31[0] + M33[7] * M31[1] + M33[8] * M31[2]);
	}

	/**
	* @brief       Matrix Norm2
	* @param[in]   r           int      The number of rows in the matrix M
	* @param[in]   c           int      The number of cols in the matrix M
	* @param[in]   M           double   Matrix
	* @return      double      Returns the two norm of a matrix
	* @note
	* @par History:
	*              
	* @internals
	*/
	inline double MatrixNorm2(int r, int c, const double M[])
	{
		int i;
		double val = 0;

		int n = (r == 1 ? c : (c == 1 ? r : 0));

		for (i = 0; i < n; i++)
		{
			val += M[i] * M[i];
		}

		return sqrt(val);
	}

	/** @brief  dot product of two vectors */
	inline double MatrixDot(int n, const double v1[], const double v2[])
	{
		double c = 0.0;
		while (--n >= 0)
			c += v1[n] * v2[n];
		return c;
	}

	inline void MatrixMultiply(int r1, int c1, const double M1[], int r2, int c2, const double M2[], double M3[], double scale = 1.0)
	{
		int i, j, k;
		double Sum;

		for (i = 0; i < r1; i++)
		{
			for (j = 0; j < c2; j++)
			{
				Sum = 0.0;

				for (k = 0; k < c1; k++)
				{
					Sum = Sum + *(M1 + i * c1 + k) * *(M2 + k * c2 + j);
				}

				*(M3 + i * c2 + j) = Sum * scale;
			}
		}
	}

	/**
	* @brief       计算直线的最近点和单位方向向量
	* @param[in]   X1			double[3]		起点
	* @param[in]   X2			double[3]		终点
	* @param[out]  X0			double[3]		最近点
	* @param[out]  N0			double[3]		单位方向向量
	* @param[out]  JX0_X1		double[3]		最近点对起点的偏导数
	* @param[out]  JX0_X2		double[3]		最近点对终点的偏导数
	* @param[out]  JN0_X1		double[3]		单位方向向量对起点的偏导数
	* @param[out]  JN0_X2		double[3]		单位方向向量对起点的偏导数
	* @return      bool
	* @note
	*              方向指向Z轴正方向
	* @par History:
	*
	* @internals
	*
	*/
	inline bool GetLineCPAndUnitNormByTwoPoint(const double X1[3], const double X2[3], double X0[3], double N0[3],
		double JX0_X1[9] = NULL, double JX0_X2[9] = NULL, double JN0_X1[9] = NULL, double JN0_X2[9] = NULL)
	{
		double	len = 0.0, M31[3] = { 0.0 };

		// 单位向量
		M31_M31(X2, X1, N0);
		len = MatrixNorm2(3, 1, N0);
		if (len < 1.0e-14)
			return false;
		M31Scale(1.0 / len, N0);
		if (N0[2] < 0.0)	M31Scale(-1.0, N0);

		// 最近点
		len = MatrixDot(3, X1, N0);
		M31EQU(N0, M31); M31Scale(len, M31);
		M31_M31(X1, M31, X0);

		len = MatrixDot(3, X0, N0);
		return true;
	}

	/**
	* @brief       获取直线的最小参数
	* @param[in]   X0			double[3]		最近点
	* @param[in]   Dir			double[3]		单位方向向量
	* @param[out]  R			double[9]		旋转矩阵
	* @param[out]  Alpha		double  		尺度参数
	* @return      bool
	* @note
	*
	* @par History:
	*
	* @internals
	*
	*/
	inline bool GetLineParamsByCPAndUnitNorm(const double X0[3], const double N0[3], double R[9], double *Alpha)
	{
		double	len = 0.0, M31[3] = { 0.0 }, M33[9] = { 0.0 };
		len = MatrixNorm2(3, 1, X0);

		// col[1]
		M33[0] = N0[0]; M33[3] = N0[1]; M33[6] = N0[2];

		// col[2]
		M31EQU(X0, M31); M31Scale(1.0 / len, M31);
		M33[1] = M31[0]; M33[4] = M31[1]; M33[7] = M31[2];

		// col[3]
		CrossM3(N0, X0, M31);
		M31Scale(1.0 / len, M31);
		M33[2] = M31[0]; M33[5] = M31[1]; M33[8] = M31[2];

		if (R)		M33EQU(M33, R);
		if (Alpha)	*Alpha = len;
		return true;
	}


	/**
	* @brief       获取直线的最小参数
	* @param[in]   X1			double[3]		起点
	* @param[in]   X2			double[3]		终点
	* @param[out]  R			double[9]		旋转矩阵
	* @param[out]  Alpha		double[3]		尺度参数
	* @return      bool
	* @note
	*
	* @par History:
	*
	* @internals
	*
	*/
	inline bool GetLineParamsByTwoPoint(const double X1[3], const double X2[3], double R[9], double *Alpha)
	{
		double	X0[3] = { 0.0 }, N0[3] = { 0.0 };
		GetLineCPAndUnitNormByTwoPoint(X1, X2, X0, N0);
		GetLineParamsByCPAndUnitNorm(X0, N0, R, Alpha);
		return true;
	}

	inline void convert_equation_to_closepoint(const float equation[4], float CP[3])
	{
		// 平面方程转最近点, 平面的法向量向外
		CP[0] = fabs(equation[3]) * equation[0];
		CP[1] = fabs(equation[3]) * equation[1];
		CP[2] = fabs(equation[3]) * equation[2];

		return;
	}

	inline void convert_equation_to_closepoint(const float equation[4], double CP[3])
	{
		// 平面方程转最近点, 平面的法向量向外
		CP[0] = fabs(equation[3]) * equation[0];
		CP[1] = fabs(equation[3]) * equation[1];
		CP[2] = fabs(equation[3]) * equation[2];

		return;
	}

	inline void convert_equation_to_closepoint(const double equation[4], float CP[3])
	{
		// 平面方程转最近点, 平面的法向量向外
		CP[0] = fabs(equation[3]) * equation[0];
		CP[1] = fabs(equation[3]) * equation[1];
		CP[2] = fabs(equation[3]) * equation[2];

		return;
	}

	inline void convert_equation_to_closepoint(const double equation[4], double CP[3])
	{
		// 平面方程转最近点, 平面的法向量向外
		CP[0] = fabs(equation[3]) * equation[0];
		CP[1] = fabs(equation[3]) * equation[1];
		CP[2] = fabs(equation[3]) * equation[2];

		return;
	}

	inline void convert_closepoint_to_equation(const float CP[3], float Equa[4])
	{
		float	CPNorm = sqrt(CP[0] * CP[0] + CP[1] * CP[1] + CP[2] * CP[2]);
		Equa[0] = CP[0] / CPNorm; Equa[1] = CP[1] / CPNorm; Equa[2] = CP[2] / CPNorm;
		Equa[3] = -(Equa[0] * CP[0] + Equa[1] * CP[1] + Equa[2] * CP[2]);

		return;
	}

	inline void convert_closepoint_to_equation(const float CP[3], double Equa[4])
	{
		double	CPNorm = sqrt(CP[0] * CP[0] + CP[1] * CP[1] + CP[2] * CP[2]);
		Equa[0] = CP[0] / CPNorm; Equa[1] = CP[1] / CPNorm; Equa[2] = CP[2] / CPNorm;
		Equa[3] = -(Equa[0] * CP[0] + Equa[1] * CP[1] + Equa[2] * CP[2]);

		return;
	}

	inline void convert_closepoint_to_equation(const double CP[3], float Equa[4])
	{
		double	CPNorm = sqrt(CP[0] * CP[0] + CP[1] * CP[1] + CP[2] * CP[2]);
		Equa[0] = CP[0] / CPNorm; Equa[1] = CP[1] / CPNorm; Equa[2] = CP[2] / CPNorm;
		Equa[3] = -(Equa[0] * CP[0] + Equa[1] * CP[1] + Equa[2] * CP[2]);

		return;
	}

	inline void convert_closepoint_to_equation(const double CP[3], double Equa[4])
	{
		double	CPNorm = sqrt(CP[0] * CP[0] + CP[1] * CP[1] + CP[2] * CP[2]);
		Equa[0] = CP[0] / CPNorm; Equa[1] = CP[1] / CPNorm; Equa[2] = CP[2] / CPNorm;
		Equa[3] = -(Equa[0] * CP[0] + Equa[1] * CP[1] + Equa[2] * CP[2]);

		return;
	}

	inline void transformClosedPoint(const float input[3], const double R[9], const double T[3], float output[3])
	{
		double	input_equ[4] = { 0.0 }, out_equ[4] = { 0.0 };
		double	Rw[9] = { 0.0 }, Norm_w[3] = { 0.0 }, Norm_wX[9] = { 0.0 };
		double	output_1[3] = { 0.0 }, output_2[3] = { 0.0 };
		double	dj = sqrt(input[0] * input[0] + input[1] * input[1] + input[2] * input[2]); // 原来的距离

		convert_closepoint_to_equation(input, input_equ);
		M33XM31(R, input_equ, out_equ); // 法向量

		// 第一部分
		M31EQU(out_equ, output_1);
		M31Scale(dj, output_1);

		// 第二部分
		MatrixMultiply(3, 1, out_equ, 1, 3, out_equ, Norm_wX);
		M33XM31(Norm_wX, T, output_2);

		output[0] = output_1[0] + output_2[0];
		output[1] = output_1[1] + output_2[1];
		output[2] = output_1[2] + output_2[2];
		return;
	}

	inline void transformClosedPoint(float CP[3], const double R[9], const double T[3])
	{
		float	CP1[3] = { 0.0 };
		transformClosedPoint(CP, R, T, CP1);
		CP[0] = CP1[0]; CP[1] = CP1[1]; CP[2] = CP1[2];
		return;
	}

	inline void transformClosedPoint(const double CPI[3], const double R[9], const double T[3], double CPO[3])
	{
		float	CP1[3] = { 0.0 }, CP2[3] = { 0.0 };
		CP1[0] = CPI[0]; CP1[1] = CPI[1]; CP1[2] = CPI[2];
		transformClosedPoint(CP1, R, T, CP2);
		CPO[0] = CP2[0]; CPO[1] = CP2[1]; CPO[2] = CP2[2];
		return;
	}

	inline void transformClosedPoint(double CP[3], const double R[9], const double T[3])
	{
		float	CP1[3] = { 0.0 }, CP2[3] = { 0.0 };
		CP1[0] = CP[0]; CP1[1] = CP[1]; CP1[2] = CP[2];
		transformClosedPoint(CP1, R, T, CP2);
		CP[0] = CP2[0]; CP[1] = CP2[1]; CP[2] = CP2[2];
		return;
	}

	// 直线(线段)观测值
	typedef struct tagPointObs
	{
		double	timestamp;			///< 当前帧时间
		int64_t globalId;			///< 当前点的全局ID
		int64_t	innerId;			///< 当前帧内的ID
		double	XYZ[3];				///< 坐标
		bool	Valid;				///< 有效标志位
		MapPoint* pMapPoint;		///< 全局系下的地图点

		tagPointObs()
		{
			globalId = innerId = -1;
			timestamp = 0.0;
			XYZ[0] = XYZ[1] = XYZ[2] = 0.0;
			Valid = false;
			pMapPoint = NULL;
		}
	} PointObs;

	typedef struct tagObsPair
	{
		int64_t	Id1;
		int64_t	Id2;
		lpostk::Frame*	pFrame1;
		lpostk::Frame*	pFrame2;

		tagObsPair()
		{
			Id1 = Id2 = -1;
			pFrame1 = pFrame2 = NULL;
		}

	} ObsPair;

	class MapPoint
	{
	public:
		MapPoint() = delete;
		MapPoint(const double XYZ[3])
		{
			this->XYZ[0] = this->XYZ[1] = this->XYZ[2] = 0.0;
			if (XYZ != NULL)
			{
				for (int i = 0; i < 3; i++)
					this->XYZ[i] = XYZ[i];
			}

			this->ID = nNextId++;
			this->Valid = true;
		}

		MapPoint(int64_t Id, const double XYZ[3])
		{
			this->XYZ[0] = this->XYZ[1] = this->XYZ[2] = 0.0;
			if (XYZ != NULL)
			{
				for (int i = 0; i < 3; i++)
					this->XYZ[i] = XYZ[i];
				Valid = true;
			}

			this->ID = Id;
			this->Valid = true;
		}

		virtual ~MapPoint(){}

		static long unsigned int nNextId;	///< 全局ID

		int64_t		ID;
		double		XYZ[3];
		bool		Valid;
		std::map<lpostk::Frame*, size_t> m_pFrames;	///< 每一帧中的索引
		std::list<std::pair<ceres::ResidualBlockId, ObsPair>> residuals_blkid; // ceres中的那个残差块
	};

	class MapLine
	{
	public:
		MapLine() = delete;
		MapLine(const double Alpha, const double R[9])
		{
			Eigen::Matrix3d	tmp;
			tmp << R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8];
			Eigen::Quaterniond qq = Eigen::Quaterniond(tmp);

			this->Alpha = Alpha;
			this->Rotation = liso::SO3d(qq);
			this->ID = nNextId++;
			this->Valid = true;
		}

		MapLine(int64_t Id, const double Alpha, const double R[9])
		{
			Eigen::Matrix3d	tmp;
			tmp << R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8];
			Eigen::Quaterniond qq = Eigen::Quaterniond(tmp);

			this->Alpha = Alpha;
			this->Rotation = liso::SO3d(qq);
			this->ID = Id;
			this->Valid = true;
		}

		virtual ~MapLine() {}

		static long unsigned int nNextId;	///< 全局ID

		int64_t		ID;
		double		Alpha;
		liso::SO3d	Rotation;
		bool		Valid;
		std::map<lpostk::Frame*, size_t> m_pFrames;	///< 每一帧中的索引
		std::vector<std::pair<ceres::ResidualBlockId, ObsPair>> residuals_blkid; // ceres中的那个残差块
	};

	class MapPole
	{
	public:
		MapPole() = delete;
		MapPole(const double Alpha, const double R[9], const double Radius = 0.0)
		{
			Eigen::Matrix3d	tmp;
			tmp << R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8];
			Eigen::Quaterniond qq = Eigen::Quaterniond(tmp);

			this->Alpha = Alpha;
			this->Rotation = liso::SO3d(qq);
			this->Radius = Radius;
			this->ID = nNextId++;
			this->Valid = true;
		}

		MapPole(int64_t Id, const double Alpha, const double R[9], const double Radius = 0.0)
		{
			Eigen::Matrix3d	tmp;
			tmp << R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8];
			Eigen::Quaterniond qq = Eigen::Quaterniond(tmp);

			this->Alpha = Alpha;
			this->Rotation = liso::SO3d(qq);
			this->Radius = Radius;
			this->ID = Id;
			this->Valid = true;
		}

		virtual ~MapPole() {}

		static long unsigned int nNextId;	///< 全局ID

		int64_t		ID;
		double		Alpha;
		liso::SO3d	Rotation;
		double		Radius;
		bool		Valid;
		std::map<lpostk::Frame*, size_t> m_pFrames;	///< 每一帧中的索引
		std::vector<std::pair<ceres::ResidualBlockId, ObsPair>> residuals_blkid; // ceres中的那个残差块
	};

	class MapPlane
	{
	public:
		MapPlane() = delete;
		MapPlane(const double CPW[3])
		{
			this->CPW[0] = this->CPW[1] = this->CPW[2] = 0.0;
			if (CPW != NULL)
			{
				for (int i = 0; i < 3; i++)
					this->CPW[i] = CPW[i];
			}

			this->ID = nNextId++;
			this->Valid = true;
		}

		MapPlane(int64_t Id, const double CPW[3])
		{
			this->CPW[0] = this->CPW[1] = this->CPW[2] = 0.0;
			if (CPW != NULL)
			{
				for (int i = 0; i < 3; i++)
					this->CPW[i] = CPW[i];
				Valid = true;
			}

			this->ID = Id;
			this->Valid = true;
		}

		virtual ~MapPlane() {}

		static long unsigned int nNextId;	///< 全局ID

		int64_t		ID;
		double		CPW[3];
		bool		Valid;
		std::map<lpostk::Frame*, size_t> m_pFrames;	///< 每一帧中的索引
		std::vector<std::pair<ceres::ResidualBlockId, ObsPair>> residuals_blkid; // ceres中的那个残差块
	};

	// 直线(线段)观测值
	typedef struct tagLineObs
	{
		double	timestamp;			///< 当前帧时间
		int64_t globalId;			///< 当前点的全局ID
		int64_t	innerId;			///< 当前帧内的ID
		double	P1[3];				///< 直线的起点
		double	P2[3];				///< 直线的终点
		double	P0[3];				///< 直线与XY平面的交点
		double	Normal[3];			///< 直线的方向
		bool	Valid;				///< 有效标志位
		std::map<lpostk::Frame*, size_t> m_pFrames;	///< 每一帧中的索引
		MapLine* pMapLine;			///< 全局系下的直线

		tagLineObs()
		{
			timestamp = 0.0;
			innerId = -1;
			for (int i = 0; i < 3; i++)
			{
				P1[i] = P2[i] = 0.0;
				P0[i] = Normal[i] = 0.0;
			}

			Valid = false;
			pMapLine = NULL;
		}
	} LineObs;

	// 直线关联
	typedef struct tagLineAssociation
	{
		int		line1_id;
		int		line2_id;
		double	Dst_center;
		double	Angle;				///< 轴线的夹角
		double	mse;
		double	DD;
		double	DD1;				///< 轴线最高点
		double	DD2;				///< 轴线最低点
		double	dP0[3];				///< 轴线与XY平面交点的坐标差值

		tagLineAssociation()
		{
			line1_id = line2_id = -1;
			Dst_center = Angle = 0.0;
			mse = 0.0;
			DD = DD1 = DD2 = 0.0;
			dP0[0] = dP0[1] = dP0[2] = 0.0;
		}
	} LineAssociation;

	// 圆柱观测值
	typedef struct tagPoleObs
	{
		double	timestamp;			///< 当前帧时间
		int64_t globalId;			///< 当前点的全局ID
		int64_t	innerId;			///< 当前帧内的ID
		double	P1[3];				///< 轴线的最低点
		double	P2[3];				///< 轴线的最高点
		double	P0[3];				///< 轴线与XY平面的交点
		double	Normal[3];			///< 轴线的方向
		double	R0;					///< 圆柱的半径
		bool	Valid;				///< 有效标志位
		std::map<lpostk::Frame*, size_t> m_pFrames;	///< 每一帧中的索引
		MapPole* pMapPole;			///< 全局系下的直线

		tagPoleObs()
		{
			timestamp = 0.0;
			innerId = -1;
			for (int i = 0; i < 3; i++)
			{
				P1[i] = P2[i] = 0.0;
				P0[i] = Normal[i] = 0.0;
			}

			R0 = 0.0;
			pMapPole = NULL;
			Valid = false;
		}

	} PoleObs;

	// 圆柱关联
	typedef struct tagPoleAssociation
	{
		int		pole1_id;
		int		pole2_id;
		double	Dst_center;
		double	Angle;				///< 轴线的夹角
		double	mse;
		double	DD;
		double	DD1;				///< 轴线最高点
		double	DD2;				///< 轴线最低点
		double	dP0[3];				///< 轴线与XY平面交点的坐标差值
		bool	Valid;				///< 有效标志位

		tagPoleAssociation()
		{
			pole1_id = pole2_id = -1;
			Dst_center = Angle = 0.0;
			mse = 0.0;
			DD = DD1 = DD2 = 0.0;
			dP0[0] = dP0[1] = dP0[2] = 0.0;
			Valid = false;
		}
	} PoleAssociation;

	// 圆柱观测值
	typedef struct tagPlaneObs
	{
		double	timestamp;			///< 当前帧时间
		int64_t globalId;			///< 当前点的全局ID
		int64_t	innerId;			///< 当前帧内的ID
		double	Center[3];			///< 中心点
		double	Equation[4];		///< 平面方程
		double	CP[3];				///< 平面最近点
		bool	Valid;				///< 有效标志位
		std::map<lpostk::Frame*, size_t> m_pFrames;	///< 每一帧中的索引
		MapPlane* pMapPlane;		///< 全局系下的直线

		tagPlaneObs()
		{
			timestamp = 0.0;
			innerId = -1;
			Center[0] = Center[1] = Center[2] = 0.0;
			Equation[0] = Equation[1] = Equation[2] = Equation[3] = 0.0;
			CP[0] = CP[1] = CP[2];
			pMapPlane = NULL;
			Valid = false;
		}

	} PlaneObs;

	typedef struct tagPlaneAssociation
	{
		int		plane1_id;
		int		plane2_id;
		double	Dst_center;
		double	Dst_closep;
		double	Angle;
		double	mse;
		double	dCP[3];

		tagPlaneAssociation()
		{
			plane1_id = plane2_id = -1;
			Dst_center = Dst_closep = Angle = 0.0;
			mse = 0.0;
			dCP[0] = dCP[1] = dCP[2] = 0.0;
		}

	} PlaneAssociation;

	class Frame
	{
	public:

		Frame()
		{
			timestamp = -1;

			for (int i = 0; i < 9; i++)
			{
				if (i < 3)
					tbn[i] = tln[i] = 0.0;

				Rbn[i] = Rln[i] = 0.0;
				if (i == 0 || i == 4 || i == 8)
					Rbn[i] = Rln[i] = 1.0;
			}
		}
		virtual ~Frame() { }

		bool read_points(const std::string filename)
		{
			Points.clear();
			FILE*	fp = fopen(filename.c_str(), "r");
			if (!fp)	return false;

			// 文件名去掉扩展名拿到时间
			char	FileName[MAXSIZE] = { 0 };
			getFileInfo(filename.c_str(), NULL, NULL, FileName, NULL);
			double	timestamp = std::stod(std::string(FileName));

			char	oneline[256] = { 0 };
			int		innerId = 0;
			while (!feof(fp))
			{
				if (!fgets(oneline, sizeof(oneline), fp))
					break;

				PointObs	point;
				if (sscanf(oneline, "%lld %lf %lf %lf",
					&point.globalId, point.XYZ + 0, point.XYZ + 1, point.XYZ + 2) != 4)
					continue;

				point.timestamp = timestamp;
				point.innerId = innerId++;
				Points.push_back(point);
			}

			fclose(fp);
			return true;
		}

		bool read_lines(const std::string filename)
		{
			Lines.clear();
			FILE*	fp = fopen(filename.c_str(), "r");
			if (!fp)	return false;

			// 文件名去掉扩展名拿到时间
			char	FileName[MAXSIZE] = { 0 };
			getFileInfo(filename.c_str(), NULL, NULL, FileName, NULL);
			double	timestamp = std::stod(std::string(FileName));

			char	oneline[256] = { 0 };
			int		innerId = 0;
			while (!feof(fp))
			{
				if (!fgets(oneline, sizeof(oneline), fp))
					break;

				LineObs	line; // x0 y0 z0 nx ny nz pxz pyz pz1 px2 py2 pz2
				if (sscanf(oneline, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
					&line.P0[0], &line.P0[1], &line.P0[2], &line.Normal[0], &line.Normal[1], &line.Normal[2],
					&line.P1[0], &line.P1[1], &line.P1[2], &line.P2[0], &line.P2[1], &line.P2[2]) != 12)
					continue;

				line.timestamp = timestamp;
				line.innerId = innerId++;
				Lines.push_back(line);
			}

			fclose(fp);
			return true;
		}

		bool read_poles(const std::string filename)
		{
			Poles.clear();
			FILE*	fp = fopen(filename.c_str(), "r");
			if (!fp)	return false;

			// 文件名去掉扩展名拿到时间
			char	FileName[MAXSIZE] = { 0 };
			getFileInfo(filename.c_str(), NULL, NULL, FileName, NULL);
			double	timestamp = std::stod(std::string(FileName));

			char	oneline[256] = { 0 };
			int		innerId = 0;
			while (!feof(fp))
			{
				if (!fgets(oneline, sizeof(oneline), fp))
					break;

				PoleObs	pole; // x0 y0 z0 nx ny nz r pxz pyz pz1 px2 py2 pz2
				if (sscanf(oneline, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
					&pole.P0[0], &pole.P0[1], &pole.P0[2], &pole.Normal[0], &pole.Normal[1], &pole.Normal[2], &pole.R0,
					&pole.P1[0], &pole.P1[1], &pole.P1[2], &pole.P2[0], &pole.P2[1], &pole.P2[2]) != 13)
					continue;

				pole.timestamp = timestamp;
				pole.innerId = innerId++;
				Poles.push_back(pole);
			}

			fclose(fp);
			return true;
		}

		bool read_planes(const std::string filename)
		{
			Planes.clear();
			FILE*	fp = fopen(filename.c_str(), "r");
			if (!fp)	return false;

			// 文件名去掉扩展名拿到时间
			char	FileName[MAXSIZE] = { 0 };
			getFileInfo(filename.c_str(), NULL, NULL, FileName, NULL);
			double	timestamp = std::stod(std::string(FileName));

			char	oneline[256] = { 0 };
			int		innerId = 0;
			while (!feof(fp))
			{
				if (!fgets(oneline, sizeof(oneline), fp))
					break;

				PlaneObs	plane;
				// cx cy xz rms a b c d cpx cpy cpz rms vcpx vcpy vcpz
				if (sscanf(oneline, "%lf %lf %lf %*f %lf %lf %lf %lf %lf %lf %lf %*f %*f %*f %*f",
					&plane.Center[0], &plane.Center[1], &plane.Center[2],
					&plane.Equation[0], &plane.Equation[1], &plane.Equation[2], &plane.Equation[3],
					&plane.CP[0], &plane.CP[1], &plane.CP[2]) != 10)
					continue;

				plane.timestamp = timestamp;
				plane.innerId = innerId++;
				Planes.push_back(plane);
			}

			fclose(fp);
			return true;
		}

		void showFrameInfo()
		{
			int		nPoint = 0, nLine = 0, nPole = 0, nPlane = 0;

			for (int i = 0; i < Points.size(); i++)
			{
				if (Points[i].pMapPoint == NULL)
					continue;

				if (Points[i].pMapPoint->m_pFrames.size() < 5)
					continue;

				nPoint++;
			}

			//for (std::map<int64_t, PointObs>::iterator iter = Points.begin(); iter != Points.end(); ++iter)
			//{
			//	if (iter->second.pMapPoint == NULL)
			//		continue;

			//	nPoint++;
			//}

			//for (std::map<int64_t, PointObs>::iterator iter = Points.begin(); iter != Points.end(); ++iter)
			//{
			//	if (iter->second.pMapPoint == NULL)
			//		continue;

			//	nPoint++;
			//}

			//for (std::map<int64_t, PlaneObs>::iterator iter = Planes.begin(); iter != Planes.end(); ++iter)
			//{
			//	if (iter->second.pMapPoint == NULL)
			//		continue;

			//	nPoint++;
			//}


			fprintf(stdout, "%.6lf %8d %8d %8d %8d\n", timestamp, nPoint, nLine, nPole, nPlane);

		}

		double	timestamp;
		double	Rbn[9];
		double	tbn[3];
		double	Rln[9];
		double	tln[3];
		std::vector<PointObs>	Points;
		std::vector<LineObs>	Lines;
		std::vector<PoleObs>	Poles;
		std::vector<PlaneObs>	Planes;
		bool	m_bValid;
		int64_t	m_StateIndex;

	};

	static bool CompareFrame(Frame* &frame1, Frame* &frame2)
	{
		double dt = frame1->timestamp - frame2->timestamp;
		if (dt < 0.0)   return true;
		return false;
	}

	static bool CompareFrameIdPair(std::pair<lpostk::Frame*, size_t> &x, std::pair<lpostk::Frame*, size_t> &y)
	{
		double dt = x.first->timestamp - y.first->timestamp;
		if (dt < 0.0)   return true;
		return false;
	}

	class FrameContainer
	{
	public:
		FrameContainer() {}
		virtual ~FrameContainer()
		{
			for (std::map<int64_t, MapPoint*>::iterator iter = MapPoints.begin(); iter != MapPoints.end(); ++iter)
			{
				delete iter->second; iter->second = NULL;
			}
			MapPoints.clear();

			for (std::map<int64_t, MapLine*>::iterator iter = MapLines.begin(); iter != MapLines.end(); ++iter)
			{
				delete iter->second; iter->second = NULL;
			}
			MapLines.clear();

			for (std::map<int64_t, MapPole*>::iterator iter = MapPoles.begin(); iter != MapPoles.end(); ++iter)
			{
				delete iter->second; iter->second = NULL;
			}
			MapPoles.clear();

			for (std::map<int64_t, MapPlane*>::iterator iter = MapPlanes.begin(); iter != MapPlanes.end(); ++iter)
			{
				delete iter->second; iter->second = NULL;
			}
			MapPlanes.clear();
		}

		bool open_points_folder(const std::string FolderPath)
		{
			char msg[MAXSIZE] = { 0 };
			fprintf(stdout, "Reading Lidar Point OBS File In %s ... \n", FolderPath.c_str());

			// 获取文件夹下所有文件
			std::vector<std::string>	m_FileNames;
			const int n = BrowPathFiles(FolderPath.c_str(), "*.pts.txt", m_FileNames);
			fprintf(stdout, "Found %d Lidar Point OBS Files In %s\n", n, FolderPath.c_str());
			if (n <= 0)	return false;

			int		nFrame = 0;
			std::sort(m_FileNames.begin(), m_FileNames.end());
			for (size_t i = 0; i < n; i++)
			{
				std::string	FileName = FolderPath + FILEPATHSEP + m_FileNames[i];
				double	sec = std::stod(m_FileNames[i]); // 文件名去掉扩展名拿到时间
				if (sec <= 0)	continue;

				char	FileTitle[MAXSIZE] = { 0 };
				sprintf(FileTitle, "%.6lf", sec);

				Frame	*frame = NULL;
				std::string	FrameKey = std::string(FileTitle);
				if (m_FramesMap.count(FrameKey) > 0)
				{
					// 当前帧已经存在了,直接追加进去就行了
					frame = m_FramesMap[FrameKey];
					if (false == frame->read_points(FileName))
						continue;

					frame->m_bValid = true;
					nFrame += 1;
				}
				else
				{
					frame = new Frame();
					if (false == frame->read_points(FileName))
					{
						delete frame;
						continue;
					}

					frame->m_bValid = true;
					nFrame += 1;

					// 当前帧的时间
					frame->timestamp = sec;
					m_Frames.push_back(frame);
					m_FramesMap.insert(std::pair<std::string, Frame*>(FrameKey, frame));
				}
			}

			fprintf(stdout, "Successfully Read %d Lidar Point OBS Files In %s\n", nFrame, FolderPath.c_str());
			return sortFrameContainer();

		}
		bool open_points_filelist(const std::string FileList, const std::string PrefixPath = "")
		{
			char msg[MAXSIZE] = { 0 };
			fprintf(stdout, "Reading Lidar Point OBS File In %s ... \n", FileList.c_str());

			std::vector<std::string> FileNames;
			if (!ReadFileListsFromListFile(FileList, FileNames, PrefixPath))
				return false;

			int		nFrame = 0;
			const size_t n = FileNames.size();
			for (size_t i = 0; i < n; i++)
			{
				// 文件名去掉扩展名拿到时间
				char	FileName[MAXSIZE] = { 0 };
				getFileInfo(FileNames[i].c_str(), NULL, NULL, FileName, NULL);
				double	sec = std::stod(std::string(FileName));
				if (sec <= 0)	continue;

				char	FileTitle[MAXSIZE] = { 0 };
				sprintf(FileTitle, "%.6lf", sec);

				Frame	*frame = NULL;
				std::string	FrameKey = std::string(FileTitle);
				if (m_FramesMap.count(FrameKey) > 0)
				{
					// 当前帧已经存在了,直接追加进去就行了
					frame = m_FramesMap[FrameKey];
					if (false == frame->read_points(FileNames[i]))
						continue;

					frame->m_bValid = true;
					nFrame += 1;
				}
				else
				{
					frame = new Frame();
					if (false == frame->read_points(FileNames[i]))
					{
						delete frame;
						continue;
					}

					frame->m_bValid = true;
					nFrame += 1;

					// 当前帧的时间
					frame->timestamp = sec;
					m_Frames.push_back(frame);
					m_FramesMap.insert(std::pair<std::string, Frame*>(FrameKey, frame));
				}
			}

			printf("Successfully Read %d Lidar Point OBS Files From %s\n", nFrame, FileList.c_str());
			return sortFrameContainer();
		}

		bool open_lines_folder(const std::string FolderPath)
		{
			char msg[MAXSIZE] = { 0 };
			fprintf(stdout, "Reading Lidar Line OBS File In %s ... \n", FolderPath.c_str());

			// 获取文件夹下所有文件
			std::vector<std::string>	m_FileNames;
			const int n = BrowPathFiles(FolderPath.c_str(), "*.lines.txt", m_FileNames);
			fprintf(stdout, "Found %d Lidar Line OBS Files In %s\n", n, FolderPath.c_str());
			if (n <= 0)	return false;

			int		nFrame = 0;
			std::sort(m_FileNames.begin(), m_FileNames.end());
			for (size_t i = 0; i < n; i++)
			{
				std::string	FileName = FolderPath + FILEPATHSEP + m_FileNames[i];
				double	sec = std::stod(m_FileNames[i]); // 文件名去掉扩展名拿到时间
				if (sec <= 0)	continue;

				char	FileTitle[MAXSIZE] = { 0 };
				sprintf(FileTitle, "%.6lf", sec);

				Frame	*frame = NULL;
				std::string	FrameKey = std::string(FileTitle);
				if (m_FramesMap.count(FrameKey) > 0)
				{
					// 当前帧已经存在了,直接追加进去就行了
					frame = m_FramesMap[FrameKey];
					if (false == frame->read_lines(FileName))
						continue;

					frame->m_bValid = true;
					nFrame += 1;
				}
				else
				{
					frame = new Frame();
					if (false == frame->read_lines(FileName))
					{
						delete frame;
						continue;
					}

					frame->m_bValid = true;
					nFrame += 1;

					// 当前帧的时间
					frame->timestamp = sec;
					m_Frames.push_back(frame);
					m_FramesMap.insert(std::pair<std::string, Frame*>(FrameKey, frame));
				}
			}

			fprintf(stdout, "Successfully Read %d Lidar Line OBS Files In %s\n", nFrame, FolderPath.c_str());
			return sortFrameContainer();
		}

		bool open_lines_filelist(const std::string FileList, const std::string PrefixPath = "")
		{
			char msg[MAXSIZE] = { 0 };
			fprintf(stdout, "Reading Lidar Line OBS File In %s ... \n", FileList.c_str());

			std::vector<std::string> FileNames;
			if (!ReadFileListsFromListFile(FileList, FileNames, PrefixPath))
				return false;

			int		nFrame = 0;
			const size_t n = FileNames.size();
			for (size_t i = 0; i < n; i++)
			{
				// 文件名去掉扩展名拿到时间
				char	FileName[MAXSIZE] = { 0 };
				getFileInfo(FileNames[i].c_str(), NULL, NULL, FileName, NULL);
				double	sec = std::stod(std::string(FileName));
				if (sec <= 0)	continue;

				char	FileTitle[MAXSIZE] = { 0 };
				sprintf(FileTitle, "%.6lf", sec);

				Frame	*frame = NULL;
				std::string	FrameKey = std::string(FileTitle);
				if (m_FramesMap.count(FrameKey) > 0)
				{
					// 当前帧已经存在了,直接追加进去就行了
					frame = m_FramesMap[FrameKey];
					if (false == frame->read_lines(FileNames[i]))
						continue;

					frame->m_bValid = true;
					nFrame += 1;
				}
				else
				{
					frame = new Frame();
					if (false == frame->read_lines(FileNames[i]))
					{
						delete frame;
						continue;
					}

					frame->m_bValid = true;
					nFrame += 1;

					// 当前帧的时间
					frame->timestamp = sec;
					m_Frames.push_back(frame);
					m_FramesMap.insert(std::pair<std::string, Frame*>(FrameKey, frame));
				}
			}

			printf("Successfully Read %d Lidar Line OBS Files From %s\n", nFrame, FileList.c_str());
			return sortFrameContainer();
		}

		bool open_poles_folder(const std::string FolderPath)
		{
			char msg[MAXSIZE] = { 0 };
			fprintf(stdout, "Reading Lidar Pole OBS File In %s ... \n", FolderPath.c_str());

			// 获取文件夹下所有文件
			std::vector<std::string>	m_FileNames;
			const int n = BrowPathFiles(FolderPath.c_str(), "*.poles.txt", m_FileNames);
			fprintf(stdout, "Found %d Lidar Pole OBS Files In %s\n", n, FolderPath.c_str());
			if (n <= 0)	return false;

			int		nFrame = 0;
			std::sort(m_FileNames.begin(), m_FileNames.end());
			for (size_t i = 0; i < n; i++)
			{
				std::string	FileName = FolderPath + FILEPATHSEP + m_FileNames[i];
				double	sec = std::stod(m_FileNames[i]); // 文件名去掉扩展名拿到时间
				if (sec <= 0)	continue;

				char	FileTitle[MAXSIZE] = { 0 };
				sprintf(FileTitle, "%.6lf", sec);

				Frame	*frame = NULL;
				std::string	FrameKey = std::string(FileTitle);
				if (m_FramesMap.count(FrameKey) > 0)
				{
					// 当前帧已经存在了,直接追加进去就行了
					frame = m_FramesMap[FrameKey];
					if (false == frame->read_poles(FileName))
						continue;

					frame->m_bValid = true;
					nFrame += 1;
				}
				else
				{
					frame = new Frame();
					if (false == frame->read_poles(FileName))
					{
						delete frame;
						continue;
					}

					frame->m_bValid = true;
					nFrame += 1;

					// 当前帧的时间
					frame->timestamp = sec;
					m_Frames.push_back(frame);
					m_FramesMap.insert(std::pair<std::string, Frame*>(FrameKey, frame));
				}
			}

			fprintf(stdout, "Successfully Read %d Lidar Pole OBS Files In %s\n", nFrame, FolderPath.c_str());
			return sortFrameContainer();
		}

		bool open_poles_filelist(const std::string FileList, const std::string PrefixPath = "")
		{
			char msg[MAXSIZE] = { 0 };
			fprintf(stdout, "Reading Lidar Pole OBS File In %s ... \n", FileList.c_str());

			std::vector<std::string> FileNames;
			if (!ReadFileListsFromListFile(FileList, FileNames, PrefixPath))
				return false;

			int		nFrame = 0;
			const size_t n = FileNames.size();
			for (size_t i = 0; i < n; i++)
			{
				// 文件名去掉扩展名拿到时间
				char	FileName[MAXSIZE] = { 0 };
				getFileInfo(FileNames[i].c_str(), NULL, NULL, FileName, NULL);
				double	sec = std::stod(std::string(FileName));
				if (sec <= 0)	continue;

				char	FileTitle[MAXSIZE] = { 0 };
				sprintf(FileTitle, "%.6lf", sec);

				Frame	*frame = NULL;
				std::string	FrameKey = std::string(FileTitle);
				if (m_FramesMap.count(FrameKey) > 0)
				{
					// 当前帧已经存在了,直接追加进去就行了
					frame = m_FramesMap[FrameKey];
					if (false == frame->read_poles(FileNames[i]))
						continue;

					frame->m_bValid = true;
					nFrame += 1;
				}
				else
				{
					frame = new Frame();
					if (false == frame->read_poles(FileNames[i]))
					{
						delete frame;
						continue;
					}

					frame->m_bValid = true;
					nFrame += 1;

					// 当前帧的时间
					frame->timestamp = sec;
					m_Frames.push_back(frame);
					m_FramesMap.insert(std::pair<std::string, Frame*>(FrameKey, frame));
				}
			}

			printf("Successfully Read %d Lidar Pole OBS Files From %s\n", nFrame, FileList.c_str());
			return sortFrameContainer();
		}

		bool open_planes_folder(const std::string FolderPath)
		{
			char msg[MAXSIZE] = { 0 };
			fprintf(stdout, "Reading Lidar Line OBS File In %s ... \n", FolderPath.c_str());

			// 获取文件夹下所有文件
			std::vector<std::string>	m_FileNames;
			const int n = BrowPathFiles(FolderPath.c_str(), "*.planes.txt", m_FileNames);
			fprintf(stdout, "Found %d Lidar Line OBS Files In %s\n", n, FolderPath.c_str());
			if (n <= 0)	return false;

			int		nFrame = 0;
			std::sort(m_FileNames.begin(), m_FileNames.end());
			for (size_t i = 0; i < n; i++)
			{
				std::string	FileName = FolderPath + FILEPATHSEP + m_FileNames[i];
				double	sec = std::stod(m_FileNames[i]); // 文件名去掉扩展名拿到时间
				if (sec <= 0)	continue;

				char	FileTitle[MAXSIZE] = { 0 };
				sprintf(FileTitle, "%.6lf", sec);

				Frame	*frame = NULL;
				std::string	FrameKey = std::string(FileTitle);
				if (m_FramesMap.count(FrameKey) > 0)
				{
					// 当前帧已经存在了,直接追加进去就行了
					frame = m_FramesMap[FrameKey];
					if (false == frame->read_planes(FileName))
						continue;

					frame->m_bValid = true;
					nFrame += 1;
				}
				else
				{
					frame = new Frame();
					if (false == frame->read_planes(FileName))
					{
						delete frame;
						continue;
					}

					frame->m_bValid = true;
					nFrame += 1;

					// 当前帧的时间
					frame->timestamp = sec;
					m_Frames.push_back(frame);
					m_FramesMap.insert(std::pair<std::string, Frame*>(FrameKey, frame));
				}
			}

			fprintf(stdout, "Successfully Read %d Lidar Plane OBS Files In %s\n", nFrame, FolderPath.c_str());
			return sortFrameContainer();

		}

		bool open_planes_filelist(const std::string FileList, const std::string PrefixPath = "")
		{
			char msg[MAXSIZE] = { 0 };
			fprintf(stdout, "Reading Lidar Plane OBS File In %s ... \n", FileList.c_str());

			std::vector<std::string> FileNames;
			if (!ReadFileListsFromListFile(FileList, FileNames, PrefixPath))
				return false;

			int		nFrame = 0;
			const size_t n = FileNames.size();
			for (size_t i = 0; i < n; i++)
			{
				// 文件名去掉扩展名拿到时间
				char	FileName[MAXSIZE] = { 0 };
				getFileInfo(FileNames[i].c_str(), NULL, NULL, FileName, NULL);
				double	sec = std::stod(std::string(FileName));
				if (sec <= 0)	continue;

				char	FileTitle[MAXSIZE] = { 0 };
				sprintf(FileTitle, "%.6lf", sec);

				Frame	*frame = NULL;
				std::string	FrameKey = std::string(FileTitle);
				if (m_FramesMap.count(FrameKey) > 0)
				{
					// 当前帧已经存在了,直接追加进去就行了
					frame = m_FramesMap[FrameKey];
					if (false == frame->read_planes(FileNames[i]))
						continue;

					frame->m_bValid = true;
					nFrame += 1;
				}
				else
				{
					frame = new Frame();
					if (false == frame->read_planes(FileNames[i]))
					{
						delete frame;
						continue;
					}

					frame->m_bValid = true;
					nFrame += 1;

					// 当前帧的时间
					frame->timestamp = sec;
					m_Frames.push_back(frame);
					m_FramesMap.insert(std::pair<std::string, Frame*>(FrameKey, frame));
				}
			}

			printf("Successfully Read %d Lidar Plane OBS Files From %s\n", nFrame, FileList.c_str());
			return sortFrameContainer();
		}

		bool sortFrameContainer()
		{
			std::list<Frame*>::iterator iter;

			// 按照时间重新排序
			m_Frames.sort(CompareFrame);

			m_nFrames = (unsigned long)m_Frames.size();
			m_bValid = (m_nFrames > 0);
			if (m_bValid == false)	return m_bValid;

			iter = m_Frames.begin();
			m_BeginTime = (*iter)->timestamp;
			if (m_BeginTime < 0)
			{
				iter++;
				m_BeginTime = (*iter)->timestamp;
			}
			m_EndTime = m_Frames.back()->timestamp;

			// 计算名义上的时间间隔
			double	last_time = -1, dt = 0.0, min_dt = 9999.999;

			iter = m_Frames.begin();
			last_time = (*iter)->timestamp; iter++;
			for (; iter != m_Frames.end(); iter++)
			{
				if (last_time >= 0)
				{
					dt = fabs((*iter)->timestamp - last_time);
					if (dt < min_dt)
						min_dt = dt;
				}

				last_time = (*iter)->timestamp;
			}
			m_nominalDeltT = min_dt;

			// 头部插入一个空的数据体,方便后续迭代器读取
			if (m_Frames.front()->timestamp != -1)
				m_Frames.push_front(new Frame());

			return m_bValid;
		}

		void showFramesInfo()
		{
			for (std::list<Frame*>::iterator iter = m_Frames.begin(); iter != m_Frames.end(); iter++)
				(*iter)->showFrameInfo();
		}

		void printMapPoints()
		{
			for (std::map<int64_t, MapPoint*>::iterator iter = MapPoints.begin(); iter != MapPoints.end(); iter++)
			{
				MapPoint*	pMapPoint = iter->second;
				fprintf(stdout, "%10lld %14.3lf %14.3lf %14.3lf\n", pMapPoint->ID, pMapPoint->XYZ[0], pMapPoint->XYZ[1], pMapPoint->XYZ[2]);
			}

			return;
		}

		void initMapPointsByAvg()
		{
			// 遍历所有地图点
			for (std::map<int64_t, MapPoint*>::iterator iter = MapPoints.begin(); iter != MapPoints.end(); ++iter)
			{
				MapPoint* pPointPoint = iter->second;

				// 观察到的所有帧
				int		count = 0;
				double	XYZ[3] = { 0.0 }, AvgXYZ[3] = { 0.0 };
				std::map<lpostk::Frame*, size_t>::iterator iter2;
				for (iter2 = iter->second->m_pFrames.begin(); iter2 != iter->second->m_pFrames.end(); iter2++)
				{
					const Frame*	pFrame = iter2->first;
					const size_t	innerId = iter2->second;
					const PointObs &point = pFrame->Points[innerId];
					if (pPointPoint->ID != point.globalId)
					{
						int hint = 1;
					}

					XYZ[0] = pFrame->Rln[0] * point.XYZ[0] + pFrame->Rln[1] * point.XYZ[1] + pFrame->Rln[2] * point.XYZ[2] + pFrame->tln[0];
					XYZ[1] = pFrame->Rln[3] * point.XYZ[0] + pFrame->Rln[4] * point.XYZ[1] + pFrame->Rln[5] * point.XYZ[2] + pFrame->tln[1];
					XYZ[2] = pFrame->Rln[6] * point.XYZ[0] + pFrame->Rln[7] * point.XYZ[1] + pFrame->Rln[8] * point.XYZ[2] + pFrame->tln[2];

					count++;
					AvgXYZ[0] += XYZ[0]; AvgXYZ[1] += XYZ[1]; AvgXYZ[2] += XYZ[2];
				}

				if (count > 0)
				{
					AvgXYZ[0] /= count; AvgXYZ[1] /= count; AvgXYZ[2] /= count;
				}

				pPointPoint->XYZ[0] = AvgXYZ[0];
				pPointPoint->XYZ[1] = AvgXYZ[1];
				pPointPoint->XYZ[2] = AvgXYZ[2];
			}

		}

		void TrackPoint()
		{
			// 遍历所有帧
			for (std::list<Frame*>::iterator iter = m_Frames.begin(); iter != m_Frames.end(); iter++)
			{
				Frame*	frame = *iter;
				if (!frame->m_bValid || frame->timestamp < 0 || frame->Points.size() <= 0)
					continue;

				// 遍历所有点
				for (size_t i = 0; i < frame->Points.size(); i++)
				{
					int64_t key = frame->Points[i].globalId;
					if (key < 0)	continue;
					if (MapPoints.count(key) > 0)
					{
						// 插入
						MapPoints[key]->m_pFrames.insert(std::pair<Frame*, size_t>(frame, i));
						frame->Points[i].pMapPoint = MapPoints[key];
					}
					else
					{
						// 新建
						double	XYZ[3] = { 0.0 };
						XYZ[0] = frame->Rln[0] * frame->Points[i].XYZ[0] + frame->Rln[1] * frame->Points[i].XYZ[1] + frame->Rln[2] * frame->Points[i].XYZ[2] + frame->tln[0];
						XYZ[1] = frame->Rln[3] * frame->Points[i].XYZ[0] + frame->Rln[4] * frame->Points[i].XYZ[1] + frame->Rln[5] * frame->Points[i].XYZ[2] + frame->tln[1];
						XYZ[2] = frame->Rln[6] * frame->Points[i].XYZ[0] + frame->Rln[7] * frame->Points[i].XYZ[1] + frame->Rln[8] * frame->Points[i].XYZ[2] + frame->tln[2];
						//fprintf(stdout, "%10lld %14.3lf %14.3lf %14.3lf\n", key, XYZ[0], XYZ[1], XYZ[2]);

						MapPoint* pMapPoint = new MapPoint(key, XYZ);
						pMapPoint->m_pFrames.insert(std::pair<Frame*, size_t>(frame, i));
						frame->Points[i].pMapPoint = pMapPoint;

						// 加入
						MapPoints.insert(std::pair<int64_t, MapPoint*>(key, pMapPoint));
					}
				}
			}

			// 初始化地图点坐标
			//initMapPointsByAvg();

			// 排序
			for (std::map<int64_t, MapPoint*>::iterator iter = MapPoints.begin(); iter != MapPoints.end(); ++iter)
			{
				MapPoint* pPointPoint = iter->second;
				fprintf(stdout, "%010lld %8zd\n", pPointPoint->ID, pPointPoint->m_pFrames.size());

				// 排序
				std::vector<std::pair<lpostk::Frame*, size_t>> vec(pPointPoint->m_pFrames.begin(), pPointPoint->m_pFrames.end());
				std::sort(vec.begin(), vec.end(), CompareFrameIdPair);

				// 重新赋值
				pPointPoint->m_pFrames.clear();
				for (size_t i = 0; i < vec.size(); i++)
					pPointPoint->m_pFrames.insert(vec[i]);
			}

			for (std::map<int64_t, MapPoint*>::iterator iter = MapPoints.begin(); iter != MapPoints.end(); ++iter)
			{
				MapPoint* pPointPoint = iter->second;
				fprintf(stdout, "%010lld %8zd\n", pPointPoint->ID, pPointPoint->m_pFrames.size());

				// 观察到的所有帧
				for (std::map<lpostk::Frame*, size_t>::iterator iter2 = iter->second->m_pFrames.begin();
					iter2 != iter->second->m_pFrames.end(); iter2++)
				{
					const size_t	innerId = iter2->second;
					const PointObs	&point = iter2->first->Points[innerId];
					if (pPointPoint->ID != point.globalId)
					{
						int hint = 1;
					}

				}
			}
		}

		void TrackLine()
		{
			// 遍历所有帧
			for (std::list<Frame*>::iterator iter = m_Frames.begin(); iter != m_Frames.end(); iter++)
			{
				Frame*	frame = *iter;
				if (!frame->m_bValid || frame->timestamp < 0 || frame->Lines.size() <= 0)
					continue;

				// 遍历所有点
				for (size_t i = 0; i < frame->Lines.size(); i++)
				{
					int64_t key = frame->Lines[i].globalId;
					if (key < 0)	continue;
					if (MapLines.count(key) > 0)
					{
						// 插入
						MapLines[key]->m_pFrames.insert(std::pair<Frame*, size_t>(frame, i));
						frame->Lines[i].pMapLine = MapLines[key];
					}
					else
					{
						// 新建
						double	X1_L[3] = { 0.0 }, X2_L[3] = { 0.0 }, X1_W[3] = { 0.0 }, X2_W[3] = { 0.0 };
						double	Alpha = 0.0, R[9] = { 0.0 };
						X1_W[0] = frame->Rln[0] * frame->Lines[i].P1[0] + frame->Rln[1] * frame->Lines[i].P1[1] + frame->Rln[2] * frame->Lines[i].P1[2] + frame->tln[0];
						X1_W[1] = frame->Rln[3] * frame->Lines[i].P1[0] + frame->Rln[4] * frame->Lines[i].P1[1] + frame->Rln[5] * frame->Lines[i].P1[2] + frame->tln[1];
						X1_W[2] = frame->Rln[6] * frame->Lines[i].P1[0] + frame->Rln[7] * frame->Lines[i].P1[1] + frame->Rln[8] * frame->Lines[i].P1[2] + frame->tln[2];
						X2_W[0] = frame->Rln[0] * frame->Lines[i].P2[0] + frame->Rln[1] * frame->Lines[i].P2[1] + frame->Rln[2] * frame->Lines[i].P2[2] + frame->tln[0];
						X2_W[1] = frame->Rln[3] * frame->Lines[i].P2[0] + frame->Rln[4] * frame->Lines[i].P2[1] + frame->Rln[5] * frame->Lines[i].P2[2] + frame->tln[1];
						X2_W[2] = frame->Rln[6] * frame->Lines[i].P2[0] + frame->Rln[7] * frame->Lines[i].P2[1] + frame->Rln[8] * frame->Lines[i].P2[2] + frame->tln[2];
						GetLineParamsByTwoPoint(X1_W, X2_W, R, &Alpha);

						MapLine* pMapLine = new MapLine(key, Alpha, R);
						pMapLine->m_pFrames.insert(std::pair<Frame*, size_t>(frame, i));
						frame->Lines[i].pMapLine = pMapLine;

						// 加入
						MapLines.insert(std::pair<int64_t, MapLine*>(key, pMapLine));
					}
				}
			}

			// 排序
			for (std::map<int64_t, MapLine*>::iterator iter = MapLines.begin(); iter != MapLines.end(); ++iter)
			{
				MapLine* pMapLine = iter->second;
				fprintf(stdout, "%010lld %8zd\n", pMapLine->ID, pMapLine->m_pFrames.size());

				// 排序
				std::vector<std::pair<lpostk::Frame*, size_t>> vec(pMapLine->m_pFrames.begin(), pMapLine->m_pFrames.end());
				std::sort(vec.begin(), vec.end(), CompareFrameIdPair);

				// 重新赋值
				pMapLine->m_pFrames.clear();
				for (size_t i = 0; i < vec.size(); i++)
					pMapLine->m_pFrames.insert(vec[i]);
			}

			for (std::map<int64_t, MapLine*>::iterator iter = MapLines.begin(); iter != MapLines.end(); ++iter)
			{
				MapLine* pMapLine = iter->second;
				fprintf(stdout, "%010lld %8zd\n", pMapLine->ID, pMapLine->m_pFrames.size());

				// 观察到的所有帧
				for (std::map<lpostk::Frame*, size_t>::iterator iter2 = iter->second->m_pFrames.begin();
					iter2 != iter->second->m_pFrames.end(); iter2++)
				{
					const size_t	innerId = iter2->second;
					const LineObs	&line = iter2->first->Lines[innerId];
					if (pMapLine->ID != line.globalId)
					{
						int hint = 1;
					}
				}
			}
		}

		void TrackPole()
		{
			// 遍历所有帧
			for (std::list<Frame*>::iterator iter = m_Frames.begin(); iter != m_Frames.end(); iter++)
			{
				Frame*	frame = *iter;
				if (!frame->m_bValid || frame->timestamp < 0 || frame->Poles.size() <= 0)
					continue;

				// 遍历所有点
				for (size_t i = 0; i < frame->Poles.size(); i++)
				{
					int64_t key = frame->Poles[i].globalId;
					if (key < 0)	continue;
					if (MapPoles.count(key) > 0)
					{
						// 插入
						MapPoles[key]->m_pFrames.insert(std::pair<Frame*, size_t>(frame, i));
						frame->Poles[i].pMapPole = MapPoles[key];
					}
					else
					{
						// 新建
						double	X1_L[3] = { 0.0 }, X2_L[3] = { 0.0 }, X1_W[3] = { 0.0 }, X2_W[3] = { 0.0 };
						double	Alpha = 0.0, R[9] = { 0.0 };
						X1_W[0] = frame->Rln[0] * frame->Poles[i].P1[0] + frame->Rln[1] * frame->Poles[i].P1[1] + frame->Rln[2] * frame->Poles[i].P1[2] + frame->tln[0];
						X1_W[1] = frame->Rln[3] * frame->Poles[i].P1[0] + frame->Rln[4] * frame->Poles[i].P1[1] + frame->Rln[5] * frame->Poles[i].P1[2] + frame->tln[1];
						X1_W[2] = frame->Rln[6] * frame->Poles[i].P1[0] + frame->Rln[7] * frame->Poles[i].P1[1] + frame->Rln[8] * frame->Poles[i].P1[2] + frame->tln[2];
						X2_W[0] = frame->Rln[0] * frame->Poles[i].P2[0] + frame->Rln[1] * frame->Poles[i].P2[1] + frame->Rln[2] * frame->Poles[i].P2[2] + frame->tln[0];
						X2_W[1] = frame->Rln[3] * frame->Poles[i].P2[0] + frame->Rln[4] * frame->Poles[i].P2[1] + frame->Rln[5] * frame->Poles[i].P2[2] + frame->tln[1];
						X2_W[2] = frame->Rln[6] * frame->Poles[i].P2[0] + frame->Rln[7] * frame->Poles[i].P2[1] + frame->Rln[8] * frame->Poles[i].P2[2] + frame->tln[2];
						GetLineParamsByTwoPoint(X1_W, X2_W, R, &Alpha);

						MapPole* pMapPole = new MapPole(key, Alpha, R);
						pMapPole->m_pFrames.insert(std::pair<Frame*, size_t>(frame, i));
						frame->Poles[i].pMapPole = pMapPole;

						// 加入
						MapPoles.insert(std::pair<int64_t, MapPole*>(key, pMapPole));
					}
				}
			}

			// 排序
			for (std::map<int64_t, MapPole*>::iterator iter = MapPoles.begin(); iter != MapPoles.end(); ++iter)
			{
				MapPole* pMapPole = iter->second;
				fprintf(stdout, "%010lld %8zd\n", pMapPole->ID, pMapPole->m_pFrames.size());

				// 排序
				std::vector<std::pair<lpostk::Frame*, size_t>> vec(pMapPole->m_pFrames.begin(), pMapPole->m_pFrames.end());
				std::sort(vec.begin(), vec.end(), CompareFrameIdPair);

				// 重新赋值
				pMapPole->m_pFrames.clear();
				for (size_t i = 0; i < vec.size(); i++)
					pMapPole->m_pFrames.insert(vec[i]);
			}

			for (std::map<int64_t, MapPole*>::iterator iter = MapPoles.begin(); iter != MapPoles.end(); ++iter)
			{
				MapPole* pMapPole = iter->second;
				fprintf(stdout, "%010lld %8zd\n", pMapPole->ID, pMapPole->m_pFrames.size());

				// 观察到的所有帧
				for (std::map<lpostk::Frame*, size_t>::iterator iter2 = iter->second->m_pFrames.begin();
					iter2 != iter->second->m_pFrames.end(); iter2++)
				{
					const size_t	innerId = iter2->second;
					const LineObs	&pole = iter2->first->Lines[innerId];
					if (pMapPole->ID != pole.globalId)
					{
						int hint = 1;
					}
				}
			}
		}

		void TrackPlane()
		{
			// 遍历所有帧
			for (std::list<Frame*>::iterator iter = m_Frames.begin(); iter != m_Frames.end(); iter++)
			{
				Frame*	frame = *iter;
				if (!frame->m_bValid || frame->timestamp < 0 || frame->Planes.size() <= 0)
					continue;

				// 遍历所有点
				for (size_t i = 0; i < frame->Planes.size(); i++)
				{
					int64_t key = frame->Planes[i].globalId;
					if (key < 0)	continue;
					if (MapPlanes.count(key) > 0)
					{
						// 插入
						MapPlanes[key]->m_pFrames.insert(std::pair<Frame*, size_t>(frame, i));
						frame->Planes[i].pMapPlane = MapPlanes[key];
					}
					else
					{
						// 新建
						double	CPW[3] = { 0.0 };
						transformClosedPoint(frame->Planes[i].CP, frame->Rln, frame->tln, CPW);

						MapPlane* pPointPlane = new MapPlane(key, CPW);
						pPointPlane->m_pFrames.insert(std::pair<Frame*, size_t>(frame, i));
						frame->Planes[i].pMapPlane = pPointPlane;

						// 加入
						MapPlanes.insert(std::pair<int64_t, MapPlane*>(key, pPointPlane));
					}
				}
			}

			// 排序
			for (std::map<int64_t, MapPlane*>::iterator iter = MapPlanes.begin(); iter != MapPlanes.end(); ++iter)
			{
				MapPlane* pMapPlane = iter->second;
				fprintf(stdout, "%010lld %8zd\n", pMapPlane->ID, pMapPlane->m_pFrames.size());

				// 排序
				std::vector<std::pair<lpostk::Frame*, size_t>> vec(pMapPlane->m_pFrames.begin(), pMapPlane->m_pFrames.end());
				std::sort(vec.begin(), vec.end(), CompareFrameIdPair);

				// 重新赋值
				pMapPlane->m_pFrames.clear();
				for (size_t i = 0; i < vec.size(); i++)
					pMapPlane->m_pFrames.insert(vec[i]);
			}

			for (std::map<int64_t, MapPlane*>::iterator iter = MapPlanes.begin(); iter != MapPlanes.end(); ++iter)
			{
				MapPlane* pMapPlane = iter->second;
				fprintf(stdout, "%010lld %8zd\n", pMapPlane->ID, pMapPlane->m_pFrames.size());

				// 观察到的所有帧
				for (std::map<lpostk::Frame*, size_t>::iterator iter2 = iter->second->m_pFrames.begin();
					iter2 != iter->second->m_pFrames.end(); iter2++)
				{
					const size_t	innerId = iter2->second;
					const PlaneObs &point = iter2->first->Planes[innerId];
					if (pMapPlane->ID != point.globalId)
					{
						int hint = 1;
					}

				}
			}
		}

		void set_frame_pose(std::shared_ptr<liso::Trajectory> &trajectory, const Eigen::Matrix3d &Rlb, const Eigen::Vector3d &tlb)
		{
			std::list<Frame*>::iterator iter;
			for (iter = m_Frames.begin(); iter != m_Frames.end(); iter++)
			{
				liso::SE3d pose = trajectory->pose((*iter)->timestamp);
				Eigen::Matrix3d Rbe = pose.unit_quaternion().toRotationMatrix();
				Eigen::Vector3d posb = pose.translation();

				// 惯导的位姿
				(*iter)->tbn[0] = posb(0); (*iter)->tbn[1] = posb(1); (*iter)->tbn[2] = posb(2);
				(*iter)->Rbn[0] = Rbe(0, 0); (*iter)->Rbn[1] = Rbe(0, 1); (*iter)->Rbn[2] = Rbe(0, 2);
				(*iter)->Rbn[3] = Rbe(1, 0); (*iter)->Rbn[4] = Rbe(1, 1); (*iter)->Rbn[5] = Rbe(1, 2);
				(*iter)->Rbn[6] = Rbe(2, 0); (*iter)->Rbn[7] = Rbe(2, 1); (*iter)->Rbn[8] = Rbe(2, 2);

				// 雷达的位姿
				Eigen::Vector3d pos1 = pose.translation() + Rbe * tlb;
				Eigen::Matrix3d Rle = Rbe * Rlb;
				(*iter)->tln[0] = pos1(0); (*iter)->tln[1] = pos1(1); (*iter)->tln[2] = pos1(2);
				(*iter)->Rln[0] = Rle(0, 0); (*iter)->Rln[1] = Rle(0, 1); (*iter)->Rln[2] = Rle(0, 2);
				(*iter)->Rln[3] = Rle(1, 0); (*iter)->Rln[4] = Rle(1, 1); (*iter)->Rln[5] = Rle(1, 2);
				(*iter)->Rln[6] = Rle(2, 0); (*iter)->Rln[7] = Rle(2, 1); (*iter)->Rln[8] = Rle(2, 2);
			}

		}

		std::list<Frame*>	m_Frames;
		std::map<std::string, Frame*> m_FramesMap;
		double		m_BeginTime;
		double		m_EndTime;
		double		m_nominalDeltT;
		size_t		m_nFrames;
		bool		m_bValid;

		std::map<int64_t, MapPoint*> MapPoints;
		std::map<int64_t, MapLine*>	 MapLines;
		std::map<int64_t, MapPole*>	 MapPoles;
		std::map<int64_t, MapPlane*> MapPlanes;
		//std::shared_ptr<liso::Trajectory> trajectory_;
	};

	class Map
	{
	public:
		Map() {}
		virtual ~Map() {}

		std::set<PointObs*>	Points;
		std::set<LineObs*>	Lines;
		std::set<PoleObs*>	Poles;
		std::set<PlaneObs*>	Planes;
	};
}


#endif // _LIDARFFSTREAM_H_H_