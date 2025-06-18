#include "CommonFunc.h"

#include <set>
#include <algorithm>
#include <direct.h>
#include <Windows.h>

void xstrmid(const char *src, const int nPos, const int nCount, char *dst)
{
	int		i;
	const char	*str;
	char	c;

	str = src + nPos;

	for (i = 0; i < nCount; i++)
	{
		c = *(str + i);
		if (c)
		{
			*(dst + i) = c;
		}
		else
		{
			// 除去末尾的'\n'
			if (dst[i - 1] == '\n') dst[i - 1] = '\0';
			*(dst + i) = '\0';
			break;
		}
	}

	*(dst + nCount) = '\0';
}


void strTrim(char *str, const int len)
{
	int		ipos1 = 0, ipos2 = 0;

	for (int i = 0; i < len; i++)
	{
		if (str[i] != ' ')
		{
			ipos1 = i;
			break;
		}
	}

	for (int i = len - 1; i >= 0; i--)
	{
		if (str[i] != ' ')
		{
			ipos2 = i;
			break;
		}
	}

	int		n = 0;
	for (int i = ipos1; i <= ipos2; i++)
	{
		str[n++] = str[i];
	}
	str[n] = '\0';
}

void strTrim(std::string& src)
{
	size_t ipos = 0, npos = 0;
	bool flag = false;

	const int len = static_cast<int>(src.size());
	for (int i = 0; i < len; i++)
	{
		if (src[i] != ' ')
		{
			ipos = i;
			break;
		}
	}

	for (int i = len - 1; i >= 0; i--)
	{
		if (src[i] != ' ')
		{
			npos = i - ipos;
			break;
		}
	}

	src = src.substr(ipos, npos + 1);
}

bool getFileInfo(const char* FilePath, char *FileTitle, char *FilterExt, char *FileName, char *dir)
{
	if (FileTitle) FileTitle[0] = '\0';
	if (FilterExt) FilterExt[0] = '\0';
	if (FileName)  FileName[0] = '\0';
	if (dir)       dir[0] = '\0';

	char oneline[100] = { '\0' };
	const char *p, *q;
	int ilen = (int)strlen(FilePath);

	if (p = strrchr(FilePath, cFILEPATHSEP))	p++;
	else p = FilePath;

	if ((p) && (q = strrchr(FilePath, '.')))
	{
		if (FileTitle)
		{
			strncpy(FileTitle, p, q - p);
			FileTitle[q - p] = '\0';
		}
		if (FilterExt)
		{
			strncpy(FilterExt, q + 1, ilen - (q - FilePath));
			FilterExt[ilen - (q - FilePath)] = '\0';
		}
		if (FileName)
		{
			strncpy(FileName, p, ilen - (p - FilePath - 1));
			FileName[ilen - (p - FilePath - 1)] = '\0';
		}
	}
	else
	{
		return false;
	}

	if ((p = strrchr(FilePath, cFILEPATHSEP)) && dir)
	{
		xstrmid(FilePath, 0, (int)(p - FilePath), dir);
	}

	return true;
}

int BrowPathFiles(const char *FolderPath, const char *FilterExt, std::vector<std::string>& Files)
{
#ifdef  _WINDOWS
	char path_ext_name[1024];

	// 先拼接得到带扩展名的路径名path_ext_name字符串
	strcpy(path_ext_name, FolderPath);
	// 判断path_name是否带'/'
	int j = (int)strlen(FolderPath);
	char ls[2] = { '\0' };
	ls[0] = FolderPath[j - 1];
	if (strcmp(ls, FILEPATHSEP) != 0) strcat(path_ext_name, FILEPATHSEP);
	strcat(path_ext_name, FilterExt);

	WIN32_FIND_DATA FileData;
	HANDLE	hSearch;
	int		nCount = 0;
	BOOL	fFinished = FALSE;

	// Start searching for path_ext_name files in the current directory. 
	hSearch = FindFirstFile(path_ext_name, &FileData); //先搜索本路径下的第一个文件FileData.cFileName
	if (hSearch != INVALID_HANDLE_VALUE)
	{
		Files.push_back(FileData.cFileName);

		nCount++;									// 添加第一个文件

		while (FindNextFile(hSearch, &FileData))
		{
			Files.push_back(FileData.cFileName);
			nCount++;
		}
		if (GetLastError() != ERROR_NO_MORE_FILES)
		{
			printf("[ERROR] Couldn't find next file.\n");
		}
		// Close the search handle. 
		if (!FindClose(hSearch))
		{
			printf("[ERROR] Couldn't close search handle.\n");
		}
	}
	else// 否则没找到本路径下的第一个文件,即没找到任何文件
	{
		printf("[ERROR] No file found.\n");
	}

	if (nCount > 1)	sort(Files.begin(), Files.end());
	return nCount;

#else
	string directory(FolderPath);
	string ext(FilterExt);

	int nCount = 0;
	int iret = 1;

	regmatch_t pm[1];
	regex_t reg;

	//生成文件扩展名对应的正则表达式

	int PosDot = ext.find(".");
	ext.insert(PosDot, "\\");

	int PosFirstStar = ext.find("*");
	int PosLastStar = ext.rfind("*");

	if (PosFirstStar == PosLastStar)
	{
		ext = ext.insert(PosFirstStar, ".");
	}
	else
	{
		ext = ext.insert(PosFirstStar, ".");
		ext = ext.insert(PosLastStar, ".");
	}

	iret = regcomp(&reg, ext.c_str(), REG_EXTENDED | REG_NEWLINE);

	if (iret != 0)
	{
		return -1;
	}

	//打开目录
	DIR *dir = opendir(FolderPath);
	if (dir == NULL)
	{
		printf("[ERROR] %s is not a directory or not exist\n", FolderPath);
		return 0;
	}

	struct dirent* d_ent = NULL;

	// linux每个目录下面都有一个"."和".."要把这两个都去掉
	char dot[3] = ".";
	char dotdot[6] = "..";

	// 一行一行的读目录下的文件,文件的属性放到dirent的变量中
	while ((d_ent = readdir(dir)) != NULL)
	{
		// 忽略 "." 和 ".."
		if ((strcmp(d_ent->d_name, dot) != 0) && (strcmp(d_ent->d_name, dotdot) != 0))
		{
			// d_type可以看到当前的东西的类型,DT_DIR代表当前都到的是目录,在usr/include/dirent.h中定义的
			if (d_ent->d_type != DT_DIR)
			{
				string d_name(d_ent->d_name);
				//printf("%s\n",d_ent->d_name);

				iret = regexec(&reg, d_name.c_str(), 1, pm, 0);

				if (iret == 0)
				{
					Files.push_back(string(d_ent->d_name));
					nCount++;
				}

			}

		}

	}

	closedir(dir);
	return nCount;
#endif

}

bool ReadFileListsFromListFile(const std::string ListFile, std::vector<std::string> &FileLists, const std::string PrefixPath)
{
	// step 1 : 获取文件列表
	FileLists.clear();
	FILE*	fp = fopen(ListFile.c_str(), "r");
	if (!fp)	return false;

	char	oneline[MAXSIZE] = { 0 }, FileName[256] = { 0 }, *p = NULL;
	double	timestamp = 0.0;
	while (!feof(fp))
	{
		if (fgets(oneline, sizeof(oneline), fp) == NULL)
			break;

		// 文件路径处理
		p = strchr(oneline, '\n');
		if (p)	*p = '\0';
		p = oneline;
		while (*p)
		{
			if (*p == '\\' || *p == '/')
				*p = cFILEPATHSEP;

			p++;
		}

		// 文件名检查
		getFileInfo(oneline, NULL, NULL, FileName, NULL);

		if (sscanf(FileName, "%lf", &timestamp) != 1)
			continue;
		//if (strstr(FileName, "pcd") == NULL && strstr(FileName, "PCD") == NULL)
		//	continue;

		std::string rawPCDFile = std::string(oneline);
		if (PrefixPath != "")
		{
			// 如果指定了前缀,去掉“./”或者".\\"
			if (strlen(oneline) > 2 && (oneline[0] == '.' && (oneline[1] == '\\' || oneline[1] == '/')))
			{
				oneline[0] = oneline[1] = ' ';
				strTrim(oneline, (int)strlen(oneline));
			}
			rawPCDFile = PrefixPath + FILEPATHSEP + std::string(oneline);
		}

		FileLists.push_back(rawPCDFile);
	}
	fclose(fp);

	// 去重复
	std::set<std::string> FileSet(FileLists.begin(), FileLists.end());
	FileLists.assign(FileSet.begin(), FileSet.end());

	// 按时间排序
	std::sort(FileLists.begin(), FileLists.end());
	return FileLists.size() > 0;
}


void XYZ2LLH(const double XYZ[3], double LLH[3], int CoorSys)
{
	double a = gs_WGS84_a;
	double e2 = gs_WGS84_e2;

	if (CoorSys == 1)
	{
		a = gs_CGCS2000_a;
		e2 = gs_CGCS2000_e2;
	}

	double X = XYZ[0];
	double Y = XYZ[1];
	double Z = XYZ[2];

	double r2 = X * X + Y * Y;
	double z = 0.0;
	double zk = 0.0;
	double v = a;
	double sinp = 0.0;

	for (z = Z, zk = 0.0; fabs(z - zk) >= 1E-4;)
	{
		zk = z;
		sinp = z / sqrt(r2 + z * z);
		v = a / sqrt(1.0 - e2 * sinp*sinp);
		z = Z + v * e2*sinp;
	}
	LLH[0] = r2 > 1E-12 ? atan(z / sqrt(r2)) : (Z > 0.0 ? M_PI / 2.0 : -M_PI / 2.0);
	LLH[1] = r2 > 1E-12 ? atan2(Y, X) : 0.0;
	LLH[2] = sqrt(r2 + z * z) - v;
}

Eigen::Vector3d XYZ2LLH(const Eigen::Vector3d XYZ, int CoorSys)
{
	double a = gs_WGS84_a;
	double e2 = gs_WGS84_e2;

	if (CoorSys == 1)
	{
		a = gs_CGCS2000_a;
		e2 = gs_CGCS2000_e2;
	}

	double X = XYZ(0);
	double Y = XYZ(1);
	double Z = XYZ(2);

	double r2 = X * X + Y * Y;
	double z = 0.0;
	double zk = 0.0;
	double v = a;
	double sinp = 0.0;

	for (z = Z, zk = 0.0; fabs(z - zk) >= 1E-4;)
	{
		zk = z;
		sinp = z / sqrt(r2 + z * z);
		v = a / sqrt(1.0 - e2 * sinp*sinp);
		z = Z + v * e2*sinp;
	}

	Eigen::Vector3d LLH;
	LLH(0) = r2 > 1E-12 ? atan(z / sqrt(r2)) : (Z > 0.0 ? M_PI / 2.0 : -M_PI / 2.0);
	LLH(1) = r2 > 1E-12 ? atan2(Y, X) : 0.0;
	LLH(2) = sqrt(r2 + z * z) - v;
	return LLH;
}


void XYZ2LLH(const Eigen::Vector3d XYZ, double LLH[3], int CoorSys)
{
	double a = gs_WGS84_a;
	double e2 = gs_WGS84_e2;

	if (CoorSys == 1)
	{
		a = gs_CGCS2000_a;
		e2 = gs_CGCS2000_e2;
	}

	double X = XYZ(0);
	double Y = XYZ(1);
	double Z = XYZ(2);

	double r2 = X * X + Y * Y;
	double z = 0.0;
	double zk = 0.0;
	double v = a;
	double sinp = 0.0;

	for (z = Z, zk = 0.0; fabs(z - zk) >= 1E-4;)
	{
		zk = z;
		sinp = z / sqrt(r2 + z * z);
		v = a / sqrt(1.0 - e2 * sinp*sinp);
		z = Z + v * e2*sinp;
	}

	LLH[0] = r2 > 1E-12 ? atan(z / sqrt(r2)) : (Z > 0.0 ? M_PI / 2.0 : -M_PI / 2.0);
	LLH[1] = r2 > 1E-12 ? atan2(Y, X) : 0.0;
	LLH[2] = sqrt(r2 + z * z) - v;
	return ;
}

void LLH2XYZ(const double LLH[3], double XYZ[3], int CoorSys)
{
	double a = gs_WGS84_a;
	double e2 = gs_WGS84_e2;

	if (CoorSys == 1)
	{
		a = gs_CGCS2000_a;
		e2 = gs_CGCS2000_e2;
	}

	double sinp = sin(LLH[0]), cosp = cos(LLH[0]), sinl = sin(LLH[1]), cosl = cos(LLH[1]);
	double v = a / sqrt(1.0 - e2 * sinp*sinp);

	XYZ[0] = (v + LLH[2])*cosp*cosl;
	XYZ[1] = (v + LLH[2])*cosp*sinl;
	XYZ[2] = (v*(1.0 - e2) + LLH[2])*sinp;
}

Eigen::Vector3d LLH2XYZ(const Eigen::Vector3d LLH, int CoorSys)
{
	double a = gs_WGS84_a;
	double e2 = gs_WGS84_e2;

	if (CoorSys == 1)
	{
		a = gs_CGCS2000_a;
		e2 = gs_CGCS2000_e2;
	}

	double sinp = sin(LLH[0]), cosp = cos(LLH[0]), sinl = sin(LLH[1]), cosl = cos(LLH[1]);
	double v = a / sqrt(1.0 - e2 * sinp*sinp);

	Eigen::Vector3d XYZ;
	XYZ[0] = (v + LLH[2])*cosp*cosl;
	XYZ[1] = (v + LLH[2])*cosp*sinl;
	XYZ[2] = (v*(1.0 - e2) + LLH[2])*sinp;
	return XYZ;
}

double getGravityLocal(const double LLH[3])
{
	double m_ra[6];
	m_ra[0] = 9.7803267715;
	m_ra[1] = 0.0052790414;
	m_ra[2] = 0.0000232718;
	m_ra[3] = -0.0000030876910891;
	m_ra[4] = 0.0000000043977311;
	m_ra[5] = 0.0000000000007211;

	double s = sin(LLH[0]);
	double s2 = SQR(s);
	double s4 = SQR(s2);

	double g = m_ra[0] * (1 + m_ra[1] * s2 + m_ra[2] * s4);
	//	g = gs_WGS84_Ge*(1+gs_WGS84_K*s2)/xsqrt(1-gs_WGS84_e2*s2);
	double dgh = (m_ra[3] + m_ra[4] * s2)*LLH[2] + m_ra[5] * SQR(LLH[2]);
	g += dgh;
	return g;
}

double getGravityLocal(const Eigen::Vector3d LLH)
{
	double m_ra[6];
	m_ra[0] = 9.7803267715;
	m_ra[1] = 0.0052790414;
	m_ra[2] = 0.0000232718;
	m_ra[3] = -0.0000030876910891;
	m_ra[4] = 0.0000000043977311;
	m_ra[5] = 0.0000000000007211;

	double s = sin(LLH(0));
	double s2 = SQR(s);
	double s4 = SQR(s2);

	double g = m_ra[0] * (1 + m_ra[1] * s2 + m_ra[2] * s4);
	//	g = gs_WGS84_Ge*(1+gs_WGS84_K*s2)/xsqrt(1-gs_WGS84_e2*s2);
	double dgh = (m_ra[3] + m_ra[4] * s2)*LLH(2) + m_ra[5] * SQR(LLH(2));
	g += dgh;
	return g;
}


bool getGravityECEF(const double LLH[3], double gravity[3])
{
	double g = getGravityLocal(LLH);
	g = -g;

	gravity[0] = cos(LLH[1])*cos(LLH[0])*g;
	gravity[1] = sin(LLH[1])*cos(LLH[0])*g;
	gravity[2] = sin(LLH[0])*g;
	return true;
}

Eigen::Vector3d getGravityECEF(const double LLH[3])
{
	double g = getGravityLocal(LLH);
	g = -g;

	Eigen::Vector3d gravity;
	gravity(0) = cos(LLH[1])*cos(LLH[0])*g;
	gravity(1) = sin(LLH[1])*cos(LLH[0])*g;
	gravity(2) = sin(LLH[0])*g;
	return gravity;
}
Eigen::Vector3d getGravityECEF(const Eigen::Vector3d LLH)
{
	double g = getGravityLocal(LLH);
	g = -g;

	Eigen::Vector3d gravity;
	gravity(0) = cos(LLH(1))*cos(LLH(0))*g;
	gravity(1) = sin(LLH(1))*cos(LLH(0))*g;
	gravity(2) = sin(LLH(0))*g;
	return gravity;
}

// Vn - n系(ENU)下的速度 Wnie - 地球自转引起的n旋转 Wnen - 地球表面弯曲而引起的n系旋转 RMh - 子午圈半径 RNh - 卯酉圈半径
bool getEarthPara(const double LLH[3], const double Vn[3], double Wnie[3], double Wnen[3], double* RMh, double* RNh, bool bVn)
{
	double tRMh = 0.0, tRNh = 0.0;

	double sl = sin(LLH[0]);
	double cl = cos(LLH[0]);
	double tl = sl / cl;
	double sl2 = SQR(sl);
	double sl4 = SQR(sl2);
	double W = sqrt(1 - gs_WGS84_e2 * sl2);

	tRMh = gs_WGS84_a * (1 - gs_WGS84_e2) / pow(W, 3) + LLH[2]; // 子午圈半径
	tRNh = gs_WGS84_a / W + LLH[2]; // 卯酉圈半径

	// 地球自转引起的n旋转
	if (Wnie)
	{
		Wnie[0] = 0.0;
		Wnie[1] = gs_WGS84_OMGE * cl;
		Wnie[2] = gs_WGS84_OMGE * sl;
	}

	// 地球表面弯曲而引起的n系旋转
	if (Vn && bVn && Wnen)
	{
		Wnen[0] = -Vn[1] / tRMh;
		Wnen[1] = Vn[0] / tRNh;
		Wnen[2] = Vn[0] / tRNh * tl;
	}

	if (RMh) *RMh = tRMh;
	if (RNh) *RNh = tRNh;
	return true;
}


void Azimuth2Attitude(const double azimuth[3], double attitude[3], bool bDeg)
{
	double yaw = azimuth[0];
	if (bDeg) yaw *= D2R;

	if (yaw <= M_PI)
	{
		yaw = -yaw;
	}
	else
	{
		yaw = 2 * M_PI - yaw;
	}

	attitude[0] = yaw;
	attitude[1] = azimuth[1];
	attitude[2] = azimuth[2];

	if (bDeg) attitude[0] *= R2D;
}

void Attitude2Azimuth(const double attitude[3], double azimuth[3], bool bDeg)
{
	double yaw = attitude[0];

	if (bDeg) yaw *= D2R;

	if (yaw < 0.0)
	{
		yaw = -yaw;
	}
	else
	{
		yaw = 2 * M_PI - yaw;
	}

	azimuth[0] = yaw;
	azimuth[1] = attitude[1];
	azimuth[2] = attitude[2];

	if (bDeg) azimuth[0] *= R2D;
}

void RotAngle2RotMatrix(const double a[3], double R[9])
{

	double ca = cos(a[0]), sa = sin(a[0]); // yaw   -  z
	double cr = cos(a[1]), sr = sin(a[1]); // pitch -  x
	double cb = cos(a[2]), sb = sin(a[2]); // roll  -  y

	// n -> b
	R[0] = ca * cb - sa * sb*sr;  R[1] = sa * cb + ca * sb*sr;  R[2] = -sb * cr;
	R[3] = -sa * cr;            R[4] = ca * cr;             R[5] = sr;
	R[6] = ca * sb + sa * cb*sr;  R[7] = sa * sb - ca * cb*sr;  R[8] = cb * cr;
}

bool RotMatrix2RotAngle(const double R[9], double a[3])
{
	if (fabs(R[2]) < IPS_EPSILON && fabs(R[8]) < IPS_EPSILON)
	{
		// pitch = ±90, 则输出错误
		printf("err: RotaModel is ZXY, and the pitch avoid ±90 !\n");

		a[0] = 0.0;
		a[1] = atan(R[5]);         // atan(R[1][2])
		a[2] = atan2l(R[6], R[0]); // atan2l(R[2][0], R[0][0]);
	}
	else
	{
		a[0] = atan2l(-R[3], R[4]);  // yaw   - z
		a[1] = asin(R[5]);           // pitch - x
		a[2] = atan2l(-R[2], R[8]);  // roll  - y
	}

	return true;
}

void Rxyz2enu(const double LLH[3], double R_E2L[9])
{
	double sinp = sin(LLH[0]), cosp = cos(LLH[0]), sinl = sin(LLH[1]), cosl = cos(LLH[1]);
	R_E2L[0] = -sinl;      R_E2L[1] = cosl;       R_E2L[2] = 0.0;
	R_E2L[3] = -sinp * cosl; R_E2L[4] = -sinp * sinl; R_E2L[5] = cosp;
	R_E2L[6] = cosp * cosl;  R_E2L[7] = cosp * sinl;  R_E2L[8] = sinp;
}

void Renu2xyz(const double LLH[3], double R_L2E[9])
{
	double sinp = sin(LLH[0]), cosp = cos(LLH[0]), sinl = sin(LLH[1]), cosl = cos(LLH[1]);
	R_L2E[0] = -sinl; R_L2E[1] = -sinp * cosl;  R_L2E[2] = cosp * cosl;
	R_L2E[3] = cosl;  R_L2E[4] = -sinp * sinl;  R_L2E[5] = cosp * sinl;
	R_L2E[6] = 0.0;   R_L2E[7] = cosp;        R_L2E[8] = sinp;
}

void M33XM31(const double M33[9], const double M31[3], double M31_[3])
{
	M31_[0] = M33[0] * M31[0] + M33[1] * M31[1] + M33[2] * M31[2];
	M31_[1] = M33[3] * M31[0] + M33[4] * M31[1] + M33[5] * M31[2];
	M31_[2] = M33[6] * M31[0] + M33[7] * M31[1] + M33[8] * M31[2];
}

void M33XM33(const double M33_1[9], const double M33_2[9], double M33_3[9])
{
	int i, j, k;
	double Sum;
	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {
		Sum = 0.0; for (k = 0; k < 3; k++) Sum = Sum + M33_1[i * 3 + k] * M33_2[k * 3 + j]; M33_3[i * 3 + j] = Sum;
	}
}

void M33XM33_R(const double M33_1[9], double M33_2[9])
{
	int i;
	double M33_3[9] = { 0.0 };
	M33XM33(M33_1, M33_2, M33_3);
	for (i = 0; i < 9; i++) M33_2[i] = M33_3[i];
}

void MatrixCopy(int r, int c, const double* src, double* dst)
{
	int n = r * c, i;
	for (i = 0; i < n; i++) dst[i] = src[i];
}

void MatrixTranspose(int r, int c, const double M[], double MT[])
{
	int i, j;
	for (i = 0; i < r; i++)
	{
		for (j = 0; j < c; j++)
		{
			MT[j*r + i] = M[i*c + j];
		}
	}
}

void MatrixTranspose(int r, int c, double M[])
{
	double* MT = (double*)malloc(r*c * sizeof(double));
	MatrixTranspose(r, c, M, MT);
	MatrixCopy(r, c, MT, M);
	free(MT);
}

void RotMatrixFromE2L(const double E_R[9], const double LLH[3], double L_R[9])
{
	Renu2xyz(LLH, L_R);
	M33XM33_R(E_R, L_R);
}

void Attitude2Rbe(const double attitude[3], const double LLH[3], double Rbe[9])
{
	double Rlb[9];
	RotAngle2RotMatrix(attitude, Rlb);
	Rxyz2enu(LLH, Rbe);
	M33XM33_R(Rlb, Rbe);
	MatrixTranspose(3, 3, Rbe);
}

void Rbe2Attitude(const double Rbe[9], const double LLH[3], double attitude[3])
{
	double Rlb[9];
	double Reb[9];
	MatrixTranspose(3, 3, Rbe, Reb);
	RotMatrixFromE2L(Reb, LLH, Rlb);
	RotMatrix2RotAngle(Rlb, attitude);
}
