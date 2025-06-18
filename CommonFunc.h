#ifndef _COMMONFUNC_H_H_
#define _COMMONFUNC_H_H_

#include <string>
#include <vector>
#include <Eigen/Core>

#ifndef MAXSIZE
#define MAXSIZE 4096
#endif

#ifdef _WINDOWS
#define FILEPATHSEP			"\\"
#define cFILEPATHSEP	    '\\'
#else
#define FILEPATHSEP			"/"
#define cFILEPATHSEP	    '/'
#endif

#define SQR(x)              ((x)*(x))                               ///< x*x
#define D2R                 (0.0174532925199432957692369076849)     ///< deg to rad
#define R2D                 (57.295779513082322864647721871734)     ///< rad to deg
#define IPS_EPSILON         2.2204460492503131e-016                 ///< epsilon == DBL_EPSILON

#define isEqual(a,b)        (fabs(a-b)<IPS_EPSILON ? true : false)  ///< a == b

static const double gs_WGS84_a = 6378137.0;                         ///< earth semimajor axis (WGS84) (m)
static const double gs_WGS84_b = 6356752.31425;                     ///< earth semimajor axis (WGS84) (m)
static const double gs_WGS84_FE = 1.0 / 298.257223563;              ///< earth flattening (WGS84)
static const double gs_WGS84_e2 = 2 * gs_WGS84_FE - SQR(gs_WGS84_FE);
static const double gs_CGCS2000_a = 6378137.0;                      ///< earth semimajor axis (WGS84) (m)
static const double gs_CGCS2000_b = 6356752.314140356;              ///< earth semimajor axis (WGS84) (m)
static const double gs_CGCS2000_FE = 1.0 / 298.257222101;           ///< earth flattening (WGS84)
static const double gs_CGCS2000_e2 = 2 * gs_CGCS2000_FE - SQR(gs_CGCS2000_FE);

static const double gs_WGS84_OMGE = 7.2921151467E-5;                ///< earth angular velocity (IS-GPS) (rad/s)
static const double gs_WGS84_Ge = 9.7803267714;                     ///< gravity at equator (m/s^2) 
static const double gs_WGS84_Gp = 9.8322011865;                     ///< gravity at polar   (m/s^2) 
static const double gs_WGS84_Gm = 9.7976446561;                     ///< Mean value (normal) gravity (m/s^2)
static const double gs_WGS84_Galph = 0.00193336138707;              ///< gravity formula constant

void xstrmid(const char *src, const int nPos, const int nCount, char *dst);
void strTrim(char *str, const int len);
void strTrim(std::string& src);

bool getFileInfo(const char* FilePath, char *FileTitle, char *FilterExt, char *FileName, char *dir);

int BrowPathFiles(const char *FolderPath, const char *FilterExt, std::vector<std::string>& Files);
bool ReadFileListsFromListFile(const std::string ListFile, std::vector<std::string> &FileLists, const std::string PrefixPath);

extern void XYZ2LLH(const double XYZ[3], double LLH[3], int CoorSys = 0);
extern Eigen::Vector3d XYZ2LLH(const Eigen::Vector3d XYZ, int CoorSys = 0);
extern void LLH2XYZ(const double LLH[3], double XYZ[3], int CoorSys = 0);
extern void XYZ2LLH(const Eigen::Vector3d XYZ, double LLH[3], int CoorSys = 0);
extern Eigen::Vector3d LLH2XYZ(const Eigen::Vector3d LLH, int CoorSys = 0);

extern double getGravityLocal(const double LLH[3]);
extern double getGravityLocal(const Eigen::Vector3d LLH);
extern double getGravityLocal(double LLH[3]);
extern double getGravityLocal(const Eigen::Vector3d LLH);
extern bool  getGravityECEF(const double LLH[3], double gravity[3]);
extern Eigen::Vector3d getGravityECEF(const double LLH[3]);
extern Eigen::Vector3d getGravityECEF(const Eigen::Vector3d LLH);
extern bool getEarthPara(const double LLH[3], const double Vn[3], double Wnie[3], double Wnen[3], double* RMh, double* RNh, bool bVn = false);

extern void Azimuth2Attitude(const double azimuth[3], double attitude[3], bool bDeg = false);
extern void Attitude2Azimuth(const double attitude[3], double azimuth[3], bool bDeg = false);

extern void MatrixTranspose(int r, int c, const double M[], double MT[]);
extern void MatrixTranspose(int r, int c, double M[]);

extern void M33XM31(const double M33[9], const double M31[3], double M31_[3]);
extern void M33XM33(const double M33_1[9], const double M33_2[9], double M33_3[9]);

extern bool RotMatrix2RotAngle(const double R[9], double a[3]);
extern void RotAngle2RotMatrix(const double a[3], double R[9]);

extern void Attitude2Rbe(const double attitude[3], const double LLH[3], double Rbe[9]);
extern void Rbe2Attitude(const double Rbe[9], const double LLH[3], double attitude[3]);


#endif // _COMMONFUNC_H_H_