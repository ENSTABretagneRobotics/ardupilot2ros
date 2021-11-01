// Simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS and Jeremy NICOLA.

#ifndef INTERVAL_H
#define INTERVAL_H

#ifdef _MSC_VER
// Disable some Visual Studio warnings.
#	ifndef CRT_SECURE_NO_DEPRECATE
#		define CRT_SECURE_NO_DEPRECATE
#	endif // CRT_SECURE_NO_DEPRECATE
#	ifndef _CRT_SECURE_NO_WARNINGS
#		define _CRT_SECURE_NO_WARNINGS
#	endif // _CRT_SECURE_NO_WARNINGS
//#	ifndef _CRT_NONSTDC_NO_WARNINGS
//#		define _CRT_NONSTDC_NO_WARNINGS
//#	endif // _CRT_NONSTDC_NO_WARNINGS
#endif // _MSC_VER

#include "iboolean.h"

#include <vector>
#include <iomanip>
#include <algorithm>
#include <limits> // For std::numeric_limits<double>.

#ifdef _MSC_VER
#ifndef UNREFERENCED_PARAMETER
#define UNREFERENCED_PARAMETER(P) (P)
#endif // !UNREFERENCED_PARAMETER
#endif // _MSC_VER

#ifdef __GNUC__
#undef UNREFERENCED_PARAMETER
#define UNREFERENCED_PARAMETER(P) (void)(P)
#endif // __GNUC__

#ifdef __BORLANDC__
#undef UNREFERENCED_PARAMETER
#define UNREFERENCED_PARAMETER(P) 
#endif // __BORLANDC__

// To avoid Visual Studio 2013 warning about overflow in floating-point constant arithmetic 
// each time INFINITY or NAN is used.
#if (_MSC_VER >= 1800)
#pragma warning(disable : 4056)
#endif // (_MSC_VER >= 1800)

// To avoid Visual Studio warnings that would happen for any project using interval. 
#ifdef _MSC_VER
#pragma warning(disable : 4996)
#endif // _MSC_VER

// The definition of M_PI and M_PI_2 is not in C++ standard, although it is often provided by compilers...
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // !M_PI
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif // !M_PI_2

#ifndef INFINITY
#define INFINITY HUGE_VAL
#endif // !INFINITY

#if defined(_MSC_VER) || defined(__BORLANDC__) 
// Used to define NAN (Not A Number).
#ifndef NAN
extern const unsigned long nan[2];
extern const double nan_double;
//#define NAN (std::numeric_limits<double>::quiet_NaN())
#define NAN nan_double
#define NAN_CONSTS_NEEDED
#endif // !NAN
#endif // defined(_MSC_VER) || defined(__BORLANDC__) 

// Infinity is denoted by oo.
#ifndef oo
//#define oo 1.0/0.0
//#define oo 1000000000.0
//#define oo INFINITY
#define oo (std::numeric_limits<double>::infinity())
#endif // !oo

// Declaration of the interval class.
class interval;

// Used to define NAI (Not An Interval).
#ifndef NAI
extern const interval nai;
#define NAI nai
#define NAI_CONST_NEEDED
#endif // !NAI

#ifndef DISABLE_USING_NAMESPACE_STD_INTERVAL_H
using namespace std;
#endif // !DISABLE_USING_NAMESPACE_STD_INTERVAL_H

// Include <QDataStream> and <QDebug> before this file to be able to use Qt specific features if you have Qt.
#ifdef QT_VERSION 
//class QDataStream;
//class QDebug;
#else
#define qDebug() std::cout
#endif // QT_VERSION 

// Deprecated.
typedef double reel;

//----------------------------------------------------------------------
// Useful real-valued functions
//----------------------------------------------------------------------
double Min(std::vector<double>& x);
double Max(std::vector<double>& x);
double Sign(const double x);
double Chi(const double a, const double b, const double c);
double Arccossin(const double x, const double y);
double Arg(const double x, const double y);
double Det(double ux, double uy, double vx, double vy);
double DistanceDirSegment(double mx, double my, double theta, double ax, double ay, double bx, double by);
void DistanceDirSegment(double& d, double& phi, double mx, double my, double theta, double ax, double ay, double bx, double by);
double DistanceDirSegments(double mx, double my, double theta, 
						   std::vector<double> ax, std::vector<double> ay, std::vector<double> bx, std::vector<double> by);
void DistanceDirSegments(double& d, double& phi, double mx, double my, double theta, 
						 std::vector<double> ax, std::vector<double> ay, std::vector<double> bx, std::vector<double> by);
#define DistanceDirCercle DistanceDirCircle
#define DistanceDirCercles DistanceDirCircles
double DistanceDirCircle(double mx, double my, double theta, double cx, double cy, double r);
void DistanceDirCircle(double& d, double& phi, double mx, double my, double theta, double cx, double cy, double r);
double DistanceDirCircles(double mx, double my, double theta, std::vector<double> cx, std::vector<double> cy, std::vector<double> r);
void DistanceDirCircles(double& d, double& phi, double mx, double my, double theta, 
						std::vector<double> cx, std::vector<double> cy, std::vector<double> r);
double DistanceDirSegmentsOrCircles(double mx, double my, double theta,
									std::vector<double> ax, std::vector<double> ay, std::vector<double> bx, std::vector<double> by,
									std::vector<double> cx, std::vector<double> cy, std::vector<double> r);
void DistanceDirSegmentsOrCircles(double& d, double& phi, double mx, double my, double theta,
								  std::vector<double> ax, std::vector<double> ay, std::vector<double> bx, std::vector<double> by,
								  std::vector<double> cx, std::vector<double> cy, std::vector<double> r);

// Used by q_in only.
class borne
{
public:
	double val;
	int ouverture;
	borne();
	borne(const double&, const int&);
	friend bool operator<(const borne &x, const borne &y);
};

class interval
{
public:
	double inf;
	double sup;
	bool isEmpty;

public:
	//----------------------------------------------------------------------
	// Constructors/destructors
	//----------------------------------------------------------------------
	interval();
	interval(const double&); // Pas const pour conversion double -> interval
	interval(const double&, const double&);
	interval(const interval&);
	//----------------------------------------------------------------------
	// Operators
	//----------------------------------------------------------------------
	interval& operator=(const interval&);
	friend interval operator+(const interval&, const interval&);
	friend interval operator-(const interval&);
	friend interval operator-(const interval&, const interval&);
	friend interval operator*(const double, const interval&);
	friend interval operator*(const interval&, const double);
	friend interval operator*(const interval&, const interval&);
	friend interval operator/(const interval&, const interval&);
	friend interval operator&(const interval&, const interval&);
	friend interval operator|(const interval&, const interval&);
	friend bool operator==(const interval&, const interval&);
	friend std::ostream& operator<<(std::ostream& os, const interval& a);
#ifdef QT_VERSION 
	//friend QDataStream& operator<<(QDataStream& s, const interval& i)
	//{
	//	s << i.inf << i.sup << i.isEmpty;
	//	return s;
	//}
	//friend QDataStream& operator>>(QDataStream& s, interval& i)
	//{
	//	s >> i.inf >> i.sup >> i.isEmpty;
	//	return s;
	//}
	friend QDebug operator<<(QDebug os, const interval& a)
	{
		if (a.isEmpty) os.nospace() << "EmptyInterval";
		else if (a.inf != a.sup)
		{ 
			os.nospace() << "[" << a.inf << ", " << a.sup << "] ";
		}
		else os.nospace() << a.inf;
		return os.space();
	}
#endif // QT_VERSION 
	//----------------------------------------------------------------------
	// Member functions
	//----------------------------------------------------------------------
	interval& Intersect(const interval& Y);
	bool IsEmpty(void) const { return isEmpty; }
};

//----------------------------------------------------------------------
// Interval-valued functions
//----------------------------------------------------------------------
interval Min(const interval&, const interval&);
interval Min(const interval&, const interval&, const interval&);
interval Max(const interval&, const interval&);
interval Max(const interval&, const interval&, const interval&);
interval Abs(const interval&);
interval Sign(const interval&);
interval Modulo(const interval& a, double x);
interval Sqr(const interval&);
interval Sqrt(const interval&);
interval InvSqrt(const interval&);
interval Exp(const interval&);
interval Log(const interval&);
interval Pow(const interval&, int);
interval Pow(const interval& x, int num, int den);
interval Power(const interval&, int);
interval PowRoot(const interval& x, int num, int den);
interval Cos(const interval&);
interval Sin(const interval&);
interval Tan(const interval&);
interval Atan(const interval&);
//Arg/Atan2, param order like double version?
interval Det(interval& ux, interval& uy, interval& vx, interval& vy);
interval Det(interval& ux, interval& uy, double& vx, double& vy);
interval Step(const interval&);
interval Parabole(const interval&, double, double, double);
interval Inter(const interval&, const interval&);
interval Inter(std::vector<interval> x);
interval Union(const interval&, const interval&);
interval Union(std::vector<interval> x);
interval InterMin(const interval&, const interval&, char);
interval Inflate(const interval&, double);
#define Enveloppe Envelope
interval Envelope(std::vector<double>& x);
//----------------------------------------------------------------------
// Other functions
//----------------------------------------------------------------------
double Inf(const interval&);
double Sup(const interval&);
double Center(const interval&);
double Width(const interval&);
double Volume(const interval&);
double Rad(const interval&);
double Marge(const interval&, const interval&);
#define ToReel ToReal
#define Todouble ToReal
double ToReal(const interval&);
double Rand(const interval&);
double Eloignement(const interval&, const interval&);
double AbsMax(const interval&);
bool OverLap(const interval&, const interval&);
bool Disjoint(const interval&, const interval&);
bool Subset(const interval&, const interval&);
bool Subset(const interval&, const interval&, double epsilon);
bool SubsetStrict(const interval& a, const interval& b);
iboolean In(const interval&, const interval&);
bool In(double, const interval&);
//----------------------------------------------------------------------
// Contractors
//----------------------------------------------------------------------
#define Cplus Cadd
void Cadd(interval& Z, interval& X, interval& Y, int dir = 0);
void Cadd(interval& Z, double x, interval& Y, int dir = 0);
void Cadd(interval& Z, interval& X, double y, int dir = 0);
void Cadd(double z, interval& X, interval& Y, int dir = 0);
#define Cmoins Csub
void Csub(interval& Z, interval& X, interval& Y, int dir = 0);
void Csub(interval& Z, double x, interval& Y, int dir = 0);
void Csub(interval& Z, interval& X, double y, int dir = 0);
void Csub(double z, interval& X, interval& Y, int dir = 0);
void Csub(interval& Y, interval& X, int dir = 0);
#define Cprod Cmul
void Cmul(interval& Z, interval& X, interval& Y, int dir = 0);
void Cmul(interval& Z, double x, interval& Y, int dir = 0);
void Cmul(interval& Z, interval& X, double y, int dir = 0);
void Cdiv(interval& Z, interval& X, interval& Y, int dir = 0);
#define Cegal Cequal
void Cequal(interval& Y, interval& X, int dir);
void Cequal(interval& Y, interval& X);
void Cmin(interval& a, interval& b, interval& c, int dir = 0);
void Cmin(interval& a, interval& b, interval& c, interval& d, int dir = 0);
void Cmin(interval& a, interval& b, interval& c, interval& d, interval& e, int dir = 0);
int Cmin(interval& a, std::vector<interval>& x, int dir = 0);
void Cmax(interval& a, interval& b, interval& c, int dir = 0);
void Cabs(interval& Y, interval& X, int dir = 0);
#define Csame_sign Csign
void Csign(interval& Y, interval& X);
void Csign(interval& Y, interval& X, int dir, double a = 0);
void Cchi(interval& F, interval& A, interval& B, interval& C);
void Cgeq(interval& Y, interval& X);
void Cinteger(interval&);
void Cboolean(interval&);
void Csqr(interval& Y, interval& X, int dir = 0);
void Csqrt(interval& Y, interval& X, int dir = 0);
void Cexp(interval& Y, interval& X, int dir = 0);
void Clog(interval& Y, interval& X, int dir = 0);
#define Cpower Cpow
void Cpow(interval& Y, interval& X, int n, int dir = 0);
void Ccos(interval& Y, interval& X, int dir = 0);
void Csin(interval& Y, interval& X, int dir = 0);
void Ctan(interval& Y, interval& X, int dir = 0);
void Catan(interval& Y, interval& X, int dir = 0);
void Csinc(interval& Y, interval& X, int dir = 0);
void Carg(interval& A, interval& X, interval& Y, int dir = 0); // Not fully implemented...
int CAngle(interval& X2, interval& Y2, interval& Theta, interval& X1, interval& Y1, bool StrongAngle); // Deprecated.
#define CNorm Cnorm
void Cnorm(interval& N, interval& X, interval& Y);
void Cnorm(interval& N, interval& X, interval& Y, interval& Z, int dir = 0);
#define Cdistance Cdist
void Cdist(interval& R, interval& X1, interval& Y1, interval& X2, interval& Y2);
#define CScal Cscal
void Cscal(interval& s, interval& ux, interval& uy, interval& vx, interval& vy);
void Cscal(interval& s, double& ux, double& uy, interval& vx, interval& vy);
#define CDet Cdet
void Cdet(interval& det, interval& ux, interval& uy, interval& vx, interval& vy, int dir = 0);
void Cdet(interval& det, double& ux, double& uy, interval& vx, interval& vy, int dir = 0);
void Cdet(interval& det, interval& ux, interval& uy, double& vx, double& vy, int dir = 0);
void Cstep(interval& Y, interval& X);
void Cstep(interval& Y, interval& X, int dir, double a = 0);
void Cramp(interval& Y, interval& X, int dir = 0, double a = 0);
void Cheaviside(interval& Y, interval& X, int dir = 0, double a = 0);
void Crect(interval& Z, interval& X, interval& Y, int dir = 0);
void Crect(interval& Y, interval& X, int dir = 0);
void Ctriangle(interval& Y, interval& X, int dir = 0);
void CDistanceDirLine(interval& dist, interval& mx, interval& my, interval& theta, 
					  double& ax, double& ay, double& bx, double& by);
int CDistanceDirSegment(interval& dist, interval& mx, interval& my, interval& theta, 
						double ax, double ay, double bx, double by, int dir = 0);
void CDistanceDirSegments(interval& distmin, interval& mx, interval& my, interval& theta, 
						  std::vector<double> ax, std::vector<double> ay, std::vector<double> bx, std::vector<double> by);
void CPointInLine(interval& mx, interval& my, double& ax, double& ay, double& bx, double& by);
#define CinSegment CPointInSegment
void CPointInSegment(interval& mx, interval& my, double ax, double ay, double bx, double by);
#define CinSegments CPointInSegments
void CPointInSegments(interval& mx, interval& my, std::vector<double> ax, std::vector<double> ay, std::vector<double> bx, std::vector<double> by);
#define CinCircle CPointInCircle
void CPointInCircle(interval& mx, interval& my, double cx, double cy, double r);
#define CinCircles CPointInCircles
void CPointInCircles(interval& mx, interval& my, std::vector<double> cx, std::vector<double> cy, std::vector<double> r, bool truth = true);
#define CinSegmentsOrCircles CPointInSegmentsOrCircles
void CPointInSegmentsOrCircles(interval& mx, interval& my, 
							   std::vector<double> ax, std::vector<double> ay, std::vector<double> bx, std::vector<double> by, 
							   std::vector<double> cx, std::vector<double> cy, std::vector<double> r);
void CPointOutsideSegment(interval& mx, interval& my, double& ax, double& ay, double& bx, double& by, bool outer);
void CPointOutsideSegments(interval& mx, interval& my, 
						   std::vector<double> ax, std::vector<double> ay, std::vector<double> bx, std::vector<double> by, bool outer);

void CPoseInSegment(interval& mx, interval& my, interval& phi, double& ax, double& ay, double& bx, double& by);
void CPoseInSegments(interval& mx, interval& my, interval& phi, 
					 std::vector<double> ax, std::vector<double> ay, std::vector<double> bx, std::vector<double> by);
void CPoseInCircle(interval& mx, interval& my, interval& phi, double& cx, double& cy, double& r);
void CPoseInCircles(interval& mx, interval& my, interval& phi, std::vector<double> cx, std::vector<double> cy, std::vector<double> r);
void CPoseInSegmentsOrCircles(interval& mx, interval& my, interval& malpha, 
							  std::vector<double> ax, std::vector<double> ay, std::vector<double> bx, std::vector<double> by, 
							  std::vector<double> cx, std::vector<double> cy, std::vector<double> r);
void CPoseTrans(interval& qx, interval& qy, interval& d, interval& px, interval& py, interval& theta);  //Go straight
void CPoseRotTrans(interval& qx, interval& qy, interval& beta, interval& phi, interval& d, interval& px, interval& py, interval& alpha);
void CPoseTransInWallsOrCircles(interval& px, interval& py, interval& alpha, interval& d, 
								std::vector<double> ax, std::vector<double> ay, std::vector<double> bx, std::vector<double> by, 
								std::vector<double> cx,std::vector<double> cy, std::vector<double> r);
void CPoseTransRotInWallsOrCircles(interval& px, interval& py, interval& alpha, interval& d, interval& psi, 
								   std::vector<double> ax, std::vector<double> ay, std::vector<double> bx, std::vector<double> by, 
								   std::vector<double> cx, std::vector<double> cy, std::vector<double> r);
void CPoseRotTransRotInWallsOrCircles(interval& px, interval& py, interval& alpha, interval& phi, interval& d, interval& psi, 
									  std::vector<double> ax, std::vector<double> ay, std::vector<double> bx, std::vector<double> by, 
									  std::vector<double> cx, std::vector<double> cy, std::vector<double> r);
void CPoseRotTransPointInWallsOrCircles(interval& px, interval& py, interval& alpha, interval& phi, interval& d, 
										std::vector<double> ax, std::vector<double> ay, std::vector<double> bx, std::vector<double> by, 
										std::vector<double> cx, std::vector<double> cy, std::vector<double> r);
void CPoseTransPointInWall(interval& px,interval& py, interval& alpha, interval& d, 
						   double ax, double ay, double bx, double by, bool truth = true);
void CPoseTransPointInWalls(interval& px,interval& py, interval& alpha, interval& d0, 
							std::vector<double>& ax, std::vector<double>& ay, std::vector<double>& bx, std::vector<double>& by, bool truth = true);
void CPoseTransPointInWallsOrCircles(interval& px,interval& py, interval& alpha, interval& d0, 
									 std::vector<double> ax,std::vector<double> ay,std::vector<double> bx,std::vector<double> by, 
									 std::vector<double> cx, std::vector<double> cy, std::vector<double> r, bool truth = true);
void CPoseTowardSegment(interval& mx, interval& my, interval& theta, 
						double& ax, double& ay, double& bx, double& by, bool truth = true);
#define Ccroisepas Cnocross
void Cnocross(interval& px, interval& py, interval& mx, interval& my, double& ax, double& ay, double& bx, double& by);
#define CPatteCroiseAucunSegment CLegCrossNoSegment
void CLegCrossNoSegment(interval& dist, interval& px, interval& py, interval& theta, 
						std::vector<double> ax, std::vector<double> ay, std::vector<double> bx, std::vector<double> by);
void CLegOnWalls(interval& dist, interval& px, interval& py, interval& theta, 
				 std::vector<double> ax, std::vector<double> ay, std::vector<double> bx, std::vector<double> by);
void CLegOnWallsOrCircles(interval& dist, interval& px, interval& py, interval& theta, 
						  std::vector<double> ax, std::vector<double> ay, std::vector<double> bx, std::vector<double> by, 
						  std::vector<double> cx, std::vector<double> cy, std::vector<double> r);
//------- Procedure de reduction elementaires sur les intervalles ----------
void Contract0(char, interval&, interval&, int);
//void Contract0 (char, interval&, interval&, int, int);
//void Contract0 (char, interval&, interval&, int, int n=0);
void Contract0(char, interval&, interval&, interval&, int);
//void Contract0 (char, interval&, double&, interval&, int); //Luc
void Contract0(char, interval&);
void ShowContraction(interval&, interval&, interval&, interval&);
void IntButterfly(interval& Y, interval Yo, interval dY, interval& X, interval Xo, int dir);
void Inter1(interval&, interval&, const interval&, const interval&, const interval&);
void Sucre(interval&, const interval&);
void Cnotin(interval& X, interval& Y, int dir = 0);
void C_q_in(interval& x, int q, std::vector<interval>& y);
//----------------------------------------------------------------------
// Other
//----------------------------------------------------------------------
void diffI(interval &x0, interval &x1, interval &c0, interval &c1);
// Primitive inclusion tests
iboolean TestDiskExists(const interval& X, const interval& Y, const interval& P1, const interval& P2, const interval& P3);
iboolean TestDiskForall(const interval& X, const interval& Y, const interval& P1, const interval& P2, const interval& P3);

#endif // !INTERVAL_H
