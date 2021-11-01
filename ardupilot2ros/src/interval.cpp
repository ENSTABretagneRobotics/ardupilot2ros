// Simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS and Jeremy NICOLA.

#include "interval.h"

#ifdef __GNUC__
// Disable some GCC warnings.
#if (__GNUC__ >= 9)
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#endif // (__GNUC__ >= 9)
#if (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 6)) || (__GNUC__ > 4))
#pragma GCC diagnostic push
#endif // (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 6)) || (__GNUC__ > 4))
#endif // __GNUC__

#if defined(_MSC_VER) || defined(__BORLANDC__) 
// Enable the use of isnan().
#include <float.h>
#ifndef isnan
#define isnan _isnan
#endif // !isnan
// Used to define NAN (Not A Number).
#ifdef NAN_CONSTS_NEEDED
const unsigned long nan[2] = {0xffffffff, 0xffffffff};
const double nan_double = -*(double*)nan;
#endif // NAN_CONSTS_NEEDED
#endif // defined(_MSC_VER) || defined(__BORLANDC__) 

#ifdef NAI_CONST_NEEDED
// Used to define NAI (Not An Interval).
const interval nai = interval();
#endif // NAI_CONST_NEEDED

const double PI = M_PI;
const double PI_2 = M_PI_2;
const double TROIS_PI_2 = 3 * PI_2;
const double CINQ_PI_2 = 5 * PI_2;
const double SEPT_PI_2 = 7 * PI_2;
const double DEUX_PI = 4 * PI_2;

using namespace std;

//----------------------------------------------------------------------
// Useful real-valued functions
//----------------------------------------------------------------------
void Min1(double& zmin0, double& zmin1, double a, double b, double c)
{
	if ((a <= b) && (a <= c)) { zmin0 = a; zmin1 = min(b, c); return; };
	if ((b <= a) && (b <= c)) { zmin0 = b; zmin1 = min(a, c); return; };
	if ((c <= b) && (c <= a)) { zmin0 = c; zmin1 = min(a, b); return; };
}
//----------------------------------------------------------------------
void Max1(double& zmax0, double& zmax1, double a, double b, double c)
{
	if ((a >= b) && (a >= c)) { zmax0 = a; zmax1 = max(b, c); return; };
	if ((b >= a) && (b >= c)) { zmax0 = b; zmax1 = max(a, c); return; };
	if ((c >= b) && (c >= a)) { zmax0 = c; zmax1 = max(a, b); return; };
}
//----------------------------------------------------------------------
double Min(vector<double>& x)
{
	double d = oo;
	for (unsigned int i = 0; i < x.size(); i++)
		d = min(x[i], d);
	return d;
}
//----------------------------------------------------------------------
double Max(vector<double>& x)
{
	double d = -oo;
	for (unsigned int i = 0; i < x.size(); i++)
		d = max(x[i], d);
	return d;
}
//----------------------------------------------------------------------
double Sign(const double x)
{
	if (x > 0) return (1);
	else return (0);
}
//----------------------------------------------------------------------
double Chi(const double a, const double b, const double c)
{
	if (a < 0) return (b);
	else return (c);
}
//----------------------------------------------------------------------
static double Power(double x, int n, int RndMode)
{
	int ChangeRndMode;     // for x < 0 and odd n
	double p, z;
	ChangeRndMode = ((x < 0.0) && (n % 2 == 1));
	if (ChangeRndMode) { z = -x; RndMode = -RndMode; }
	else     z = x;
	p = 1.0;
	switch (RndMode) 
	{                             // Separate while-loops used
	case -1:   
		while (n > 0)
		{
			if (n % 2 == 1) p = p*z;   //--------------------------
			n = n / 2; if (n > 0) z = z*z;
		}
		break;
	case +1:   
		while (n > 0)
		{
			if (n % 2 == 1) p = p*z;
			n = n / 2;
			if (n > 0) z = z*z;
		}
		break;
	}
	if (ChangeRndMode) return -p; else return p;
}
//----------------------------------------------------------------------
double Arccossin(const double x, const double y)
{
	if (y > 0) return (acos(x));
	else return (-acos(x));
}
//----------------------------------------------------------------------
double Arg(const double x, const double y)
{
	double r = sqrt(x*x + y*y);
	if (r == 0) return 0;
	return Arccossin(x / r, y / r);
}
//----------------------------------------------------------------------
double Det(double ux, double uy, double vx, double vy)
{
	return (ux*vy - vx*uy);
}
//----------------------------------------------------------------------
double DistanceDirSegment(double mx, double my, double theta, double ax, double ay, double bx, double by)
{      
	// Distance directionnelle du point m au segment [a,b]. La direction est donnee par theta
	double ma_x = ax - mx;
	double ma_y = ay - my;
	double mb_x = bx - mx;
	double mb_y = by - my;
	double ab_x = bx - ax;
	double ab_y = by - ay;
	double ux = cos(theta);
	double uy = sin(theta);
	double z1 = Det(ma_x, ma_y, ux, uy);
	double z2 = Det(ux, uy, mb_x, mb_y);
	double z3 = Det(ma_x, ma_y, ab_x, ab_y);
	double z4 = Det(ux, uy, ab_x, ab_y);
	double z5 = min(z1, min(z2, z3));
	double d1 = 0;
	if (z4 == 0) d1 = oo;
	else d1 = z3 / z4;
	return (Chi(z5, oo, d1));
}
//----------------------------------------------------------------------
void DistanceDirSegment(double& d,double& phi, double mx, double my, double theta, double ax, double ay, double bx, double by)
{      
	// Distance directionnelle du point m au segment [a,b].
	double ma_x=ax-mx;
	double ma_y=ay-my;
	double mb_x=bx-mx;
	double mb_y=by-my;
	double ab_x=bx-ax;
	double ab_y=by-ay;
	double ux=cos(theta);
	double uy=sin(theta);
	double z1=Det(ma_x,ma_y,ux,uy);
	double z2=Det(ux,uy,mb_x,mb_y);
	double z3=Det(ma_x,ma_y,ab_x,ab_y);
	double z4=Det(ux,uy,ab_x,ab_y);
	double z5=min(z1,min(z2,z3));
	double d1=z3/z4;
	d=Chi(z5,oo,d1);
	phi=atan2(-ab_x,ab_y); //phi is the angle of the normal vector of [a,b]
}
//----------------------------------------------------------------------
double DistanceDirSegments(double mx, double my, double theta, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{      
	// Distance directionnelle relativement a un polygone
	vector<double> dist(ax.size());
	for (unsigned int j = 0; j < ax.size(); j++)
		dist[j] = DistanceDirSegment(mx, my, theta, ax[j], ay[j], bx[j], by[j]);
	double distmin = Min(dist);
	return (distmin);
}
//----------------------------------------------------------------------
void DistanceDirSegments(double& d, double& phi, double mx, double my, double theta, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{     
	// Distance directionnelle relativement a un polygone (le triedre m-a-b doit etre direct, sinon cette distance est infinie)
	d = oo; phi = 0;
	for (unsigned int j = 0; j < ax.size(); j++)
	{
		double dj = 0, phij = 0;
		DistanceDirSegment(dj,phij,mx,my,theta,ax[j],ay[j],bx[j],by[j]);
		if (dj < d) { d=dj;phi=phij; };
	}
}
//----------------------------------------------------------------------
double DistanceDirCircle(double mx, double my, double theta, double cx, double cy, double r)
{      
	// Distance directionnelle du point m au cercle de centre c et de rayon r.
	double ux, uy, alpha, beta, a, b, c, delta, px1, py1, px2, py2, d1, d2;
	ux = cos(theta);
	uy = sin(theta);
	if (fabs(uy) > 0.00)  //pour eviter la division par zero. Il conviendrait de traiter le cas uy=0
	{
		alpha = ux / uy;
		beta = mx - my*alpha;
		a = alpha*alpha + 1;
		b = 2 * alpha*(beta - cx) - 2 * cy;
		c = (beta - cx)*(beta - cx) + cy*cy - r*r;
		delta = b*b - 4 * a*c;
		if (delta < 0) return(oo);
		py1 = (-b - sqrt(delta)) / (2 * a);
		px1 = alpha*py1 + beta;
		py2 = (-b + sqrt(delta)) / (2 * a);
		px2 = alpha*py2 + beta;
		d1 = Chi((px1 - mx)*ux + (py1 - my)*uy, oo, sqrt((px1 - mx)*(px1 - mx) + (py1 - my)*(py1 - my)));
		d2 = Chi((px2 - mx)*ux + (py2 - my)*uy, oo, sqrt((px2 - mx)*(px2 - mx) + (py2 - my)*(py2 - my)));
		return min(d1, d2);
	}
	return oo;
}
//----------------------------------------------------------------------
void DistanceDirCircle(double& d, double& phi, double mx, double my, double theta, double cx, double cy, double r)
{      
	//  d is the directional distance from m to the circle of center c and radius r.
	//  phi is the angle of the normal vector of of the impact point
	double ux,uy,alpha,beta,a,b,c,delta,px1,py1,px2,py2,d1,d2,phi1,phi2;
	ux=cos(theta);
	uy=sin(theta);
	if (fabs(uy)<0.01)    // pour eviter une division par zero, on permutte x et y
	{ double aux; aux=mx; mx=my; my=aux;  aux=cx; cx=cy; cy=aux;    aux=ux; ux=uy; uy=aux; }
	alpha=ux/uy;
	beta=mx-my*alpha;
	a=alpha*alpha+1;
	b=2*alpha*(beta-cx)-2*cy;
	c=(beta-cx)*(beta-cx)+cy*cy-r*r;
	delta=b*b-4*a*c;
	if (delta<0) {d=oo;phi=-999;return;};
	py1=(-b-sqrt(delta))/(2*a);
	px1=alpha*py1+beta;
	phi1=atan2(cy-py1,cx-px1);
	py2=(-b+sqrt(delta))/(2*a);
	px2=alpha*py2+beta;
	phi2=atan2(cy-py2,cx-px2);

	d1=oo;d2=oo;
	if ((px1-mx)*ux+(py1-my)*uy>0)   d1=sqrt((px1-mx)*(px1-mx)+(py1-my)*(py1-my));
	if ((px2-mx)*ux+(py2-my)*uy>0)   d2=sqrt((px2-mx)*(px2-mx)+(py2-my)*(py2-my));
	if (d1<d2) { d=d1; phi=phi1; } else { d=d2; phi=phi2; }
}
//----------------------------------------------------------------------
double DistanceDirCircles(double mx, double my, double theta, vector<double> cx, vector<double> cy, vector<double> r)
{      
	// Distance directionnelle relativement a plusieurs cercles.
	vector<double> dist(cx.size());
	for (unsigned int j = 0; j < cx.size(); j++)
		dist[j] = DistanceDirCircle(mx, my, theta, cx[j], cy[j], r[j]);
	double distmin = Min(dist);
	return (distmin);
}
//----------------------------------------------------------------------
void DistanceDirCircles(double& d,double& phi, double mx, double my, double theta, vector<double> cx, vector<double> cy, vector<double> r)
{       
	d = oo; phi = 0;
	for (unsigned int j = 0; j < cx.size(); j++)
	{
		double dj = 0, phij = 0;
		DistanceDirCircle(dj,phij,mx,my,theta,cx[j],cy[j],r[j]);
		if (dj < d) { d=dj; phi=phij; };
	}
}
//----------------------------------------------------------------------
double DistanceDirSegmentsOrCircles(double mx, double my, double theta,
									vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by,
									vector<double> cx, vector<double> cy, vector<double> r)
{     
	double d1a, d1b;
	d1a = DistanceDirSegments(mx, my, theta, ax, ay, bx, by);
	d1b = DistanceDirCircles(mx, my, theta, cx, cy, r);
	return min(d1a, d1b);
}
//----------------------------------------------------------------------
void DistanceDirSegmentsOrCircles(double& d, double& phi, double mx, double my, double theta,
								  vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by,
								  vector<double> cx, vector<double> cy, vector<double> r)
{     
	// returns the distance and orientation collected by a laser rangefinder in a room made with segments and circles
	double phi1a = 0, phi1b = 0, d1a = oo, d1b = oo;
	DistanceDirSegments(d1a,phi1a,mx,my,theta,ax,ay,bx,by);
	DistanceDirCircles(d1b,phi1b,mx,my,theta,cx,cy,r);
	if (d1a < d1b) { d = d1a; phi = phi1a-theta; } else { d = d1b; phi = phi1b-theta; }
}
//----------------------------------------------------------------------
borne::borne()
{
	val = 0; ouverture = 0;
}
//----------------------------------------------------------------------
borne::borne(const double& a, const int& b)
{
	val = a; ouverture = b;
}
//----------------------------------------------------------------------
bool operator<(const borne &x, const borne &y)
{
	return x.val < y.val;
}
//----------------------------------------------------------------------
// Constructors/destructors
//----------------------------------------------------------------------
interval::interval()
{
	inf = NAN; sup = NAN;
	isEmpty = true;
}
//----------------------------------------------------------------------
interval::interval(const double& m)
{
	if (m == NAN)
	{
		inf = NAN; sup = NAN;
		isEmpty = true;
	}
	else
	{
		inf = m; sup = m;
		isEmpty = false;
	}
}
//----------------------------------------------------------------------
interval::interval(const double& a, const double& b)
{
	if ((a == NAN)||(b == NAN))
	{
		inf = NAN; sup = NAN;
		isEmpty = true;
	}
	else
	{
		if (a >= b) { inf = b; sup = a; } //modif sur >=
		else { inf = a; sup = b; }
		isEmpty = false;
	}
}
//----------------------------------------------------------------------
interval::interval(const interval& a) { *this = a; }
//----------------------------------------------------------------------
// Operators
//----------------------------------------------------------------------
interval& interval::operator=(const interval& v)
{
	(*this).inf = v.inf; (*this).sup = v.sup; (*this).isEmpty = v.isEmpty;
	return *this;
}
//----------------------------------------------------------------------
interval operator+(const interval& x, const interval& y)
{
	if (x.isEmpty || y.isEmpty) return interval();
	else return (interval(x.inf + y.inf, x.sup + y.sup));
}
//----------------------------------------------------------------------
interval operator-(const interval& a)
{
	if (a.isEmpty) return interval();
	else return (interval(-(a.sup), -(a.inf)));
}
//----------------------------------------------------------------------
interval operator-(const interval& x, const interval& y)
{
	if (x.isEmpty || y.isEmpty) return interval();
	else return (x + (-y));
}
//----------------------------------------------------------------------
interval operator*(const double x, const interval& y)
{
	if (y.isEmpty) return interval();
	else return interval(x*y.inf, x*y.sup);
}
//----------------------------------------------------------------------
interval operator*(const interval& x, const double y)
{
	return interval(y*x);
}
//----------------------------------------------------------------------
interval operator*(const interval& x, const interval& y)
{
	if (x.isEmpty || y.isEmpty) return interval();
	double x1 = x.inf, y1 = y.inf, x2 = x.sup, y2 = y.sup;
	if (x1 >= 0)
	{
		if (y1 >= 0) return interval(x1*y1, x2*y2);
		else
			if (y2 <= 0) return interval(x2*y1, x1*y2);
			else return interval(x2*y1, x2*y2);
	}
	else
		if (x2 <= 0)
		{
			if (y1 >= 0) return interval(x1*y2, x2*y1);
			else if (y2 <= 0) return interval(x2*y2, x1*y1);
			else return interval(x1*y2, x1*y1);
		}
		else
		{
			if (y1 >= 0) return interval(x1*y2, x2*y2);
			else if (y2 <= 0) return interval(x2*y1, x1*y1);
			else return interval(min(x1*y2, x2*y1), max(x1*y1, x2*y2));
		}
}
//----------------------------------------------------------------------
interval operator/(const interval& a, const interval& b)
{
	if (a.isEmpty || b.isEmpty) return interval();
	if (b.inf > 0)
	{
		if (a.inf >= 0) return interval(a.inf / b.sup, a.sup / b.inf);
		else if (a.sup < 0)  return interval(a.inf / b.inf, a.sup / b.sup);
		else  return interval(a.inf / b.inf, a.sup / b.inf);
	}
	else if (b.sup < 0)
	{
		if (a.inf >= 0) return interval(a.sup / b.sup, a.inf / b.inf);
		else if (a.sup < 0)  return interval(a.sup / b.inf, a.inf / b.sup);
		else return interval(a.sup / b.sup, a.inf / b.sup);
	}
	else return interval(-oo, oo);
}
//----------------------------------------------------------------------
interval operator&(const interval& x, const interval& y)
{
	return (Inter(x, y));
}
//----------------------------------------------------------------------
interval operator|(const interval& x, const interval& y)
{
	return (Union(x, y));
}
//----------------------------------------------------------------------
bool operator==(const interval& a, const interval& b)
{
	if ((a.inf == b.inf)&&(a.sup == b.sup)&&(a.isEmpty == b.isEmpty)) return true;
	return false;
}
//----------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const interval& a)
{
	if (a.isEmpty) os << "EmptyInterval";
	else if (a.inf != a.sup)
	{
		os << "[" << setprecision(4) << a.inf << ", " << setprecision(4) << a.sup << "] ";
	}
	else os << a.inf;
	return os;
}
//----------------------------------------------------------------------
// Member functions
//----------------------------------------------------------------------
interval& interval::Intersect(const interval& Y)
{
	interval X = *this;
	interval Z = Inter(X, Y);
	*this = Z;
	return *this;
}
//----------------------------------------------------------------------
// Interval-valued functions
//----------------------------------------------------------------------
interval Min(const interval& x, const interval& y)
{
	if (x.isEmpty || y.isEmpty) return interval();
	double a = min(x.sup, y.sup);
	double b = min(x.inf, y.inf);
	if ((x.sup == x.inf) && (y.sup == y.inf)) return interval(b);
	return interval(a, b);
}
//----------------------------------------------------------------------
interval Min(const interval& x, const interval& y, const interval& z)
{
	return (Min(Min(x, y), z));
}
//----------------------------------------------------------------------
interval Max(const interval& x, const interval& y)
{
	return (-Min(-x, -y));
}
//----------------------------------------------------------------------
interval Max(const interval& x, const interval& y, const interval& z)
{
	return (Max(Max(x, y), z));
}
//----------------------------------------------------------------------
interval Abs(const interval& a)
{
	if (a.isEmpty) return interval();
	double a1 = a.inf, a2 = a.sup;
	if ((a1 >= 0) || (a2 <= 0))   return interval(fabs(a1), fabs(a2));
	if (fabs(a1)>fabs(a2))  return interval(0, fabs(a1));
	else                    return interval(0, fabs(a2));
}
//----------------------------------------------------------------------
interval Sign(const interval& x)
{
	if (x.isEmpty) return interval();
	if (x.inf > 0) return interval(1);
	if (x.sup < 0) return interval(-1);
	return interval(-1, 1);
}
//----------------------------------------------------------------------
interval Modulo(const interval& a, double x)
{
	if ((a.inf >= 0) && (a.inf<x)) return (a);
	int k = (int)floor(a.inf / x);
	double offset = x * k;
	return interval(a.inf - offset, a.sup - offset);
}
//----------------------------------------------------------------------
interval Sqr(const interval& a)
{
	double a1 = a.inf, a2 = a.sup;
	if (a.isEmpty) return interval();
	if ((a1 >= 0) || (a2 <= 0))   return interval(a1*a1, a2*a2);
	if (fabs(a1) > fabs(a2))  return interval(0, a1*a1);
	else                    return interval(0, a2*a2);
}
//----------------------------------------------------------------------
interval Sqrt(const interval& a)
{
	if ((a.isEmpty) || (a.sup < 0)) return interval();
	double a1 = a.inf, a2 = a.sup;
	if (a1 >= 0)  return (interval(sqrt(a1), sqrt(a2)));
	else return (interval(0, sqrt(a2)));
}
//----------------------------------------------------------------------
interval InvSqrt(const interval& a)
{
	if ((a.isEmpty) || (a.sup < 0)) return interval();
	double a2 = max(a.inf, a.sup);
	//if (In(a1,interval(a2))) return sqrt(a1); else
	return interval(-sqrt(a2), sqrt(a2));
}
//----------------------------------------------------------------------
interval Exp(const interval& a)
{
	if (a.isEmpty) return(interval());
	else return(interval(exp(a.inf), exp(a.sup)));
}
//----------------------------------------------------------------------
interval Log(const interval& a)
{
	if ((a.isEmpty) || (a.sup <= 0)) return interval();
	interval b = Abs(a);
	if (In(0, a)) return interval(-oo, log(a.sup));
	return interval(log(b.inf), log(b.sup));
}
//----------------------------------------------------------------------
interval Pow(const interval& x, int n)
{
	// [x]^n
	if (x.isEmpty) return interval();
	double a1 = min(x.inf, x.sup), a2 = max(x.inf, x.sup);
	if (n > 0)
	{
		if (n % 2 == 0)
		{
			if (a1*a2 < 0) return interval(0, pow(a2, n));
			else return interval(pow(a1, n), pow(a2, n));
		}
		else return interval(pow(a1, n), pow(a2, n));
	}
	else return interval(1.0) / Pow(x, abs(n));
}
//----------------------------------------------------------------------------
interval Pow(const interval& x, int num, int den)
{
	// [x]^(num/den)
	if (x.isEmpty) return interval();
	if (num*den > 0)
	{
		double a1 = min(x.inf, x.sup), a2 = max(x.inf, x.sup);
		double n = num, m = den;
		if (den % 2 == 0)
		{
			if ((a1 >= 0) || (a1*a2 <= 0))
				return interval(-pow(a2, n / m), pow(a2, n / m));
			else return interval();
		}
		else
		{
			if (a1*a2 <= 0)
				return interval(-pow(fabs(a1), n / m), pow(a2, n / m));
			if (a1 > 0) return interval(pow(a1, n / m), pow(a2, n / m));
			else return interval(-pow(fabs(a1), n / m), -pow(fabs(a2), n / m));
		}
	}
	else return interval(1.0) / Pow(x, abs(num), abs(den));
}
//----------------------------------------------------------------------
interval Power(const interval& x, int n)
{
	int  m;   double Lower, Upper;
	if (n == 0) return(interval(1.0, 1.0));
	if (n > 0)  m = n;  else  m = -n;
	if ((0.0 < x.inf) || (m % 2 == 1))
	{
		Lower = Power(x.inf, m, -1);    Upper = Power(x.sup, m, +1);
	}
	else if (0.0 > x.sup)
	{
		Lower = Power(x.sup, m, -1);    Upper = Power(x.inf, m, +1);
	}
	else  { Lower = 0.0; Upper = Power(AbsMax(x), m, +1); }
	if (n > 0) return(interval(Lower, Upper));
	else   return(1.0 / interval(Lower, Upper));    // if 0 in 'x'.
}
//----------------------------------------------------------------------------
interval PowRoot(const interval& x, int num, int den)
{
	// [x]^(num/den)
	if (x.isEmpty) return interval();
	if (num*den > 0)
	{
		double a1 = min(x.inf, x.sup), a2 = max(x.inf, x.sup), n = num, m = den;
		if (den % 2 == 0)
		{
			if (a1 >= 0) return interval(pow(a1, n / m), pow(a2, n / m));
			if (a1*a2 <= 0) return interval(0, pow(a2, n / m));
			else return interval();
		}
		else
		{
			if (a1*a2 <= 0) return interval(-pow(fabs(a1), n / m), pow(a2, n / m));
			if (a1 > 0) return interval(pow(a1, n / m), pow(a2, n / m));
			else return interval(-pow(fabs(a1), n / m), -pow(fabs(a2), n / m));
		}
	}
	else return interval(1.0) / Pow(x, abs(num), abs(den));
}
//----------------------------------------------------------------------
interval Cos(const interval& a)
{
	if (a.isEmpty) return interval();
	else return Sin(a + PI_2);
}
//----------------------------------------------------------------------
interval Sin(const interval& a)
{
	if (a.isEmpty) return interval();
	if (a.sup == a.inf) return sin(a.inf);
	interval b;
	double sin1, sin2, r1, r2;
	b = Modulo(a, DEUX_PI);
	if (Width(a) > DEUX_PI) return (interval(-1, 1));
	sin1 = sin(b.inf);   sin2 = sin(b.sup);
	if ((b.inf < TROIS_PI_2) && (b.sup > TROIS_PI_2)) r1 = -1.0;
	else if ((b.inf < SEPT_PI_2) && (b.sup > SEPT_PI_2)) r1 = -1.0;
	else r1 = ((sin1 < sin2) ? sin1 : sin2);
	if ((b.inf < PI_2) && (b.sup > PI_2)) r2 = 1.0;
	else if ((b.inf < CINQ_PI_2) && (b.sup > CINQ_PI_2)) r2 = 1.0;
	else r2 = ((sin1 > sin2) ? sin1 : sin2);
	return Inter(interval(-1, 1), interval(r1, r2));
}
//----------------------------------------------------------------------
interval Tan(const interval& a) { return Sin(a) / Cos(a); }
//----------------------------------------------------------------------
interval Atan(const interval& x)
{
	if (x.isEmpty) return interval();
	else return interval(atan(x.inf), atan(x.sup));
}
//----------------------------------------------------------------------
interval Det(interval& ux, interval& uy, interval& vx, interval& vy)
{
	return (ux*vy - vx*uy);
}
//----------------------------------------------------------------------
interval Det(interval& ux, interval& uy, double& vx, double& vy)
{
	return (vy*ux - vx*uy);
}
//----------------------------------------------------------------------
interval Step(const interval& X)
{ 
	if (X.isEmpty) return interval();
	if (X.inf > 0) return (interval(1));
	if (X.sup < 0) return (interval(0));
	return (interval(0,1));
}
//----------------------------------------------------------------------
interval Parabole(const interval& x, double a, double b, double c)
{
	return (a*Sqr(x + (b / (2 * a))) - (b*b) / (4 * a) + c);
}
//----------------------------------------------------------------------
interval Inter(const interval& a, const interval& b)
{
	if (a.isEmpty || b.isEmpty) return interval();
	interval r(0, 0);
	if ((a.inf > b.sup) || (b.inf > a.sup)) return interval();
	if (a.inf <= b.inf) r.inf = b.inf; else r.inf = a.inf;
	if (a.sup >= b.sup) r.sup = b.sup; else r.sup = a.sup;
	return r;
}
//------------------------------------------------------------------------------
interval Inter(vector<interval> x)
{
	interval r = x[0];
	for (unsigned int i = 1; i < x.size(); i++)
		r = Inter(r, x[i]);
	return r;
}
//----------------------------------------------------------------------
interval Union(const interval& a, const interval& b)
{
	if (a.isEmpty) return b;
	if (b.isEmpty) return a;
	interval r(0, 0);
	if (a.inf <= b.inf) r.inf = a.inf; else r.inf = b.inf;
	if (a.sup >= b.sup) r.sup = a.sup; else r.sup = b.sup;
	return r;
}
//----------------------------------------------------------------------
interval Union(vector<interval> x)
{
	interval r = x[0];
	for (unsigned int i = 1; i<x.size(); i++)
		r = Union(r, x[i]);
	return r;
}
//----------------------------------------------------------------------
interval InterMin(const interval& a, const interval& b, char c)
{
	interval y(0, 0);
	if (c == '-')  { interval temp(a.inf, oo); y = Inter(temp, b); }
	else { interval temp(0, a.sup); y = Inter(temp, b); }
	return y;
}
//----------------------------------------------------------------------
interval Inflate(const interval& a, double eps)
{
	interval r(a.inf - eps, a.sup + eps);
	r.isEmpty = a.isEmpty;
	return interval(r);
}
//----------------------------------------------------------------------
interval Envelope(vector<double>& x)
{
	if (x.size() == 0) return (interval());
	return (interval(Min(x), Max(x)));
}
//----------------------------------------------------------------------
// Other functions
//----------------------------------------------------------------------
double Inf(const interval& a)
{
	if (a.inf<a.sup) return a.inf;
	else return a.sup;
}
//----------------------------------------------------------------------
double Sup(const interval& a)
{
	if (a.sup>a.inf) return a.sup;
	else return a.inf;
}
//-----------------------------------------------------------------------
double Center(const interval& a)
{
	if (a.isEmpty) return oo;
	else  return (a.sup + a.inf) / 2;
}
//----------------------------------------------------------------------
double Width(const interval& a)
{
	if (a.isEmpty) return (-1);
	else return (a.sup - a.inf);
}
//----------------------------------------------------------------------
double Volume(const interval& a)
{
	if (a.isEmpty) return (-1);
	else return (a.sup - a.inf);
}
//----------------------------------------------------------------------
double Rad(const interval& a)
{
	if (a.isEmpty) return (-1);
	else return (a.sup - a.inf)/2.0;
}
//----------------------------------------------------------------------
double Marge(const interval& a, const interval& b)
{
	return min(a.inf - b.inf, b.sup - a.sup);
}
//----------------------------------------------------------------------
double ToReal(const interval& a)
{
	if ((a.isEmpty)||(a.sup != a.inf)) return NAN;
	return a.inf;
}
//----------------------------------------------------------------------
double Rand(const interval& a)
{
	if (a.isEmpty) return oo;
	else  return a.inf + (a.sup - a.inf)*rand() / RAND_MAX;
}
//----------------------------------------------------------------------
double Eloignement(const interval& x, const interval& y)
{
	if ((x.isEmpty) || (y.isEmpty)) return oo;
	double r = max(0.0, x.inf - y.sup);
	r = max(r, y.inf - x.sup);
	return r;
}
//----------------------------------------------------------------------
double AbsMax(const interval& x)              // Absolute maximum of
{
	if (x.isEmpty) return 0;
	double a = fabs(x.inf), b = fabs(x.sup);
	return max(a, b);
}
//----------------------------------------------------------------------
bool OverLap(const interval& a, const interval& b)
{
	return (a.sup>b.inf) && (b.sup > a.inf);
}
//----------------------------------------------------------------------
bool Disjoint(const interval& a, const interval& b)
{
	if (a.isEmpty || b.isEmpty) return true;
	return ((a.sup < b.inf) || (b.sup < a.inf));
}
//----------------------------------------------------------------------
bool Subset(const interval& a, const interval& b)
{
	if (a.isEmpty) return true;
	if (b.isEmpty) return false;
	return ((a.inf >= b.inf) && (a.sup <= b.sup));
}
//----------------------------------------------------------------------
bool Subset(const interval& a, const interval& b, double epsilon)
{
	if (a.isEmpty || b.isEmpty) return false;
	else return (Subset(a, b) && ((b.inf < a.inf - epsilon) || (a.sup < b.sup - epsilon)));
}
//----------------------------------------------------------------------
bool SubsetStrict(const interval& a, const interval& b)
{
	if (a.isEmpty) return true;
	if (b.isEmpty) return false;
	return ((a.inf > b.inf) && (a.sup < b.sup));
}
//----------------------------------------------------------------------
iboolean In(const interval& F, const interval& Y)
{
	interval Z = Inflate(Y, 1e-17);
	if (Disjoint(F, Z)) return false;
	if (Subset(F, Z)) return true;
	return iboolean(iperhaps);
}
//----------------------------------------------------------------------
bool In(double a, const interval& b)
{
	interval z = Inflate(b, 1e-6);
	if (b.isEmpty) return false;
	return ((z.inf <= a)&&(a <= z.sup));
}
//----------------------------------------------------------------------
// Contractors
//----------------------------------------------------------------------
void Cadd(interval& Z, interval& X, interval& Y, int dir)
{
	// Z=X+Y         =>  dir=1;
	// X=Z-Y Y=Z-X;  =>  dir=-1;
	if (dir != -1) { Z = Inter(Z, X + Y); }
	if (dir != 1) { X = Inter(X, Z - Y); Y = Inter(Y, Z - X); }
}
//----------------------------------------------------------------------
void Cadd(interval& Z, double x, interval& Y, int dir)
{
	if (dir != -1) { Z = Inter(Z, x + Y); }
	if (dir != 1) { Y = Inter(Y, Z - x); }
}
//----------------------------------------------------------------------
void Cadd(interval& Z, interval& X, double y, int dir)
{
	if (dir != -1) { Z = Inter(Z, X + y); }
	if (dir != 1) { X = Inter(X, Z - y); }
}
//----------------------------------------------------------------------
void Cadd(double z, interval& X, interval& Y, int dir)
{ 
	if (dir != 1) { X = Inter(X, z - Y); Y = Inter(Y, z - X); }
}
//----------------------------------------------------------------------
void Csub(interval& Z, interval& X, interval& Y, int dir)
{
	// Z=X-Y           =>  dir=1;
	// X=Z+Y Y=X-Z;    =>  dir=-1;
	if (dir != -1) { Z = Inter(Z, X - Y); }
	if (dir != 1) { X = Inter(X, Z + Y); Y = Inter(Y, X - Z); }
}
//----------------------------------------------------------------------
void Csub(interval& Z, double x, interval& Y, int dir)
{
	if (dir != -1) { Z = Inter(Z, x - Y); }
	if (dir != 1) { Y = Inter(Y, x - Z); }
}
//----------------------------------------------------------------------
void Csub(interval& Z, interval& X, double y, int dir)
{
	if (dir != -1) { Z = Inter(Z, X - y); }
	if (dir != 1) { X = Inter(X, Z + y); }
}
//----------------------------------------------------------------------
void Csub(double z, interval& X, interval& Y, int dir)
{
	if (dir != 1) { X = Inter(X, z + Y); Y = Inter(Y, X - z); }
}
//----------------------------------------------------------------------
void Csub(interval& Y, interval& X, int dir)
{
	// Y=-X           =>  dir=1;
	// X=-Y           =>  dir=-1; 
	if (dir != -1) Y = Inter(Y, -X);
	if (dir != 1) X = Inter(X, -Y);
}
//----------------------------------------------------------------------
void Cmul(interval& Z, interval& X, interval& Y, int dir)
{
	// Z=X*Y           =>  dir=1;
	// X=Z/Y; Y=Z/X    =>  dir=-1; 
	if (dir != -1)  
	{ 
		Z = Inter(Z, X*Y); 
	}
	if (dir != 1) 
	{
		//modifs ???

		if (In(0, Z) == false)
		{
			interval Xd, Xg;
			interval Yd, Yg;
			Xg = Inter(X, interval(-oo, 0)), Xd = Inter(X, interval(0, oo));
			Yg = Inter(Y, interval(-oo, 0)), Yd = Inter(Y, interval(0, oo));

			if (Z.inf > 0)
			{
				if (Yd.sup*Xd.sup < Z.inf)
				{
					Xd = interval();
					Yd = interval();
				}
				else
				{
					if ((Yd.sup != 0) && (Xd.sup != 0))
					{
						Yd = Inter(Yd, interval(Z.inf / Xd.sup, oo));
						Xd = Inter(Xd, interval(Z.inf / Yd.sup, oo));
					}
				}
				if (Yg.inf*Xg.inf < Z.inf)
				{
					Xg = interval();
					Yg = interval();
				}
				else
				{
					if ((Yg.inf != 0) && (Xg.inf != 0))
					{
						Yg = Inter(Yg, interval(Z.inf / Xg.inf, -oo));
						Xg = Inter(Xg, interval(Z.inf / Yg.inf, -oo));
					}
				}
				Y = Inter(Y, Union(Yd, Yg));
				X = Inter(X, Union(Xd, Xg));
			}
			if (Z.sup<0)
			{
				if (Yg.inf*Xd.sup>Z.sup)
				{
					Xd = interval();
					Yg = interval();
				}
				else
				{
					if ((Yg.inf != 0) && (Xd.sup != 0))
					{
						Yg = Inter(Yg, interval(Z.sup / Xd.sup, -oo));
						Xd = Inter(Xd, interval(Z.sup / Yg.inf, oo));
					}
				}
				if (Yd.sup*Xg.inf > Z.sup)
				{
					Yd = interval();
					Xg = interval();
				}
				else
				{
					if ((Yd.sup != 0) && (Xg.inf != 0))
					{
						Yd = Inter(Yd, interval(Z.sup / Xg.inf, oo));
						Xg = Inter(Xg, interval(Z.sup / Yd.sup, -oo));
					}
				}
				Y = Inter(Y, Union(Yd, Yg));
				X = Inter(X, Union(Xd, Xg));
			}
		}
		//modifs  ???

		X = Inter(X, Z/Y); Y = Inter(Y, Z/X); //X = Inter(X,Z/Y);
	}
}
//----------------------------------------------------------------------
/*void Cmul(interval& Z, interval& X, interval& Y, int dir)
{ 
// Z=X*Y           =>  dir=1;
// X=Z/Y; Y=Z/X    =>  dir=-1;  
if (dir != -1) { Z = Inter(Z,X*Y); }
if (dir != 1) { X = Inter(X,Z/Y); Y = Inter(Y,Z/X); X = Inter(X,Z/Y); }
}*/
//----------------------------------------------------------------------
void Cmul(interval& Z, double x, interval& Y, int dir)
{
	if (dir != -1) { Z = Inter(Z, x*Y); }
	if (dir != 1) { Y = Inter(Y, Z/x); }
}
//----------------------------------------------------------------------
void Cmul(interval& Z, interval& X, double y, int dir)
{
	if (dir != -1) { Z = Inter(Z, X*y); }
	if (dir != 1) { X = Inter(X, Z/y); }
}
//----------------------------------------------------------------------
void Cdiv(interval& Z, interval& X, interval& Y, int dir)
{
	// Z=X/Y           =>  dir=1;
	// X=Z*Y; Y=X/Z    =>  dir=-1;
	if (dir != -1) { Z = Inter(Z, X/Y); }
	if (dir != 1) { X = Inter(X, Z*Y); Y = Inter(Y, X/Z); }
}
//----------------------------------------------------------------------
void Cequal(interval& Y, interval& X, int dir)
{
	// Y=X           =>  dir=1;
	// X=Y           =>  dir=-1; 
	if (dir != -1) Y = Inter(Y, X);
	if (dir != 1) X = Inter(X, Y);
}
//----------------------------------------------------------------------
void Cequal(interval& Y, interval& X)
{
	Y = Inter(Y, X); X = Y;
}
//----------------------------------------------------------------------
void Cmin(interval& a, interval& b, interval& c, int dir)
{   
	// a=min(b,c)                       =>  dir=1;
	// b=min-1(a,c); c=min-1(a,b)       =>  dir=-1; 
	//if ((a.isEmpty||b.isEmpty)||c.isEmpty) a=b=c=interval();
	if (dir != -1) { a = Inter(a, Min(b, c)); }
	if (dir != 1) 
	{
		if (Disjoint(a, b)) c = Inter(c, a);
		else { if (Disjoint(a, c)) b = Inter(b, a); }
		interval temp(a.inf, oo);
		b = Inter(b, temp); c = Inter(c, temp);
	}
}
//----------------------------------------------------------------------
void Cmin(interval& a, interval& b, interval& c, interval& d, int dir)
{
	// contrainte quaternaire   a=min(b,c,d)
	interval z1 = Min(b, c);
	Cmin(a, z1, d, 1);
	if (dir == 1) return;
	Cmin(a, z1, d, -1);
	Cmin(z1, b, c, -1);
}
//----------------------------------------------------------------------
void Cmin(interval& a, interval& b, interval& c, interval& d, interval& e, int dir)
{   
	interval z1 = Min(b,c,d);
	Cmin(a, z1, e, 1);
	if (dir == 1) return;
	Cmin(a, z1, e, -1);
	Cmin(z1, b, c, d, -1);
}
//----------------------------------------------------------------------
int Cmin(interval& a, vector<interval>& x, int dir)
{
	vector<interval> z(x.size());
	z[0] = x[0];
	for (unsigned int i = 1; i < x.size(); i++)
	{
		z[i] = interval(-oo, oo);
		Cmin(z[i], x[i], z[i - 1], 1);
	}
	Cequal(a, z[x.size() - 1]);
	if (dir == 1) return 1;
	for (int i = (int)(x.size() - 1); i >= 1; i--)
	{
		Cmin(z[i], x[i], z[i - 1], -1);;
		if (z[i].isEmpty) return -1;
	}
	return 1;
}
//----------------------------------------------------------------------
void Cmax(interval& a, interval& b, interval& c, int dir)
{
	// a=max(b,c)                       =>  dir=1;
	// b=max-1(a,c); c=max-1(a,b)       =>  dir=-1; 
	if (dir != -1) { a = Inter(a, Max(b, c)); }
	if (dir != 1) 
	{
		if (Disjoint(a, b)) c = Inter(c, a);
		else { if (Disjoint(a, c)) b = Inter(b, a); }
		interval temp(-oo, a.sup);
		b = Inter(b, temp); c = Inter(c, temp);
	}
}
//----------------------------------------------------------------------
void Cabs(interval& Y, interval& X, int dir)
{
	// Y=|X|=max(X,-X)     =>     dir=1
	// X=Abs-1(Y)          =>     dir=-1
	interval Xd = Inter(X, interval(0, oo)), Xg = Inter(X, interval(-oo, 0));
	if (dir != -1) { interval Yd = Inter(Y, Xd), Yg = Inter(Y, -Xg); Y = Union(Yd, Yg); }
	if (dir != 1) { Xd = Inter(Xd, Y); Xg = Inter(Xg, -Y); X = Union(Xd, Xg); }
}
//----------------------------------------------------------------------
void Csign(interval& Y, interval& X)
{
	// max(A1,A2)>0
	// A1=min(X,Y)
	// A2=min(X,Y)
	interval A1(0, oo), A2(-oo, oo), A3(-oo, oo);
	interval Xmoins(0), Ymoins(0);
	Xmoins = -X; Ymoins = -Y;
	Cmin(A2, X, Y, 1);
	Cmin(A3, Xmoins, Ymoins, 1);
	Cmax(A1, A2, A3, 1);
	Cmax(A1, A2, A3, -1);
	Cmin(A3, Xmoins, Ymoins, -1);
	Y = Inter(Y, -Ymoins);
	X = Inter(X, -Xmoins);
	Cmin(A2, X, Y, -1); Cmin(A2, X, Y, 1); Cmin(A2, X, Y, -1);
}
//----------------------------------------------------------------------
void Csign(interval& Y, interval& X, int dir, double a)
{
	// sign(X)=2*step(X)-1
	interval un(1), deux(2);
	if (dir != -1)
	{
		interval A1(-oo, oo); Cstep(A1, X, 1, a);
		Y = Inter(Y, (deux*A1) - un);
		//interval Z=Inter(X,interval(0));
		//if (!Z.isEmpty) Y=Union(Y,interval(0));
	}
	if (dir != 1)
	{
		interval A2 = (Y + un) / deux;
		Cstep(A2, X, -1, a);
	}
}
//----------------------------------------------------------------------
void Cchi(interval& F, interval& A, interval& B, interval& C)
{
	if (A.sup < 0) { Cequal(B, F); }
	else if (A.inf > 0) { Cequal(C, F); };
	if (Disjoint(F, B)) { A = Inter(A, interval(0, oo)); };
	if (Disjoint(F, C)) { A = Inter(A, interval(-oo, 0)); };
	F = Union(Inter(F, B), Inter(F, C));
}
//----------------------------------------------------------------------
void Cgeq(interval& Y, interval& X)
{
	if (Y.inf >= X.sup) return;
	interval Z = Y - X;
	Z = Inter(Z, interval(0, oo));
	Csub(Z, Y, X, -1);
}
//----------------------------------------------------------------------
void Cinteger(interval& X)
{
	interval A1(ceil(X.inf), floor(X.sup));
	X = Inter(X, A1);
}
//----------------------------------------------------------------------
void Cboolean(interval& X)
{
	interval nul(0), one(1);
	interval A1 = Inter(nul, X);
	interval A2 = Inter(one, X);
	interval A = Union(A1, A2);
	X = Inter(X, A);
}
//----------------------------------------------------------------------
void Csqr(interval& Y, interval& X, int dir)
{
	// Y=sqr(X)           =>  dir=1;
	// X=sqr-1(Y)         =>  dir=-1; 
	interval Yd, Yg, Xd = Inter(X, interval(0, oo)), Xg = Inter(X, interval(-oo, 0));
	if (!Xd.isEmpty) { Xd = Inter(Xd, Sqrt(Y));  Yd = Inter(Y, Sqr(Xd)); }
	if (!Xg.isEmpty) { Xg = Inter(Xg, -Sqrt(Y)); Yg = Inter(Y, Sqr(Xg)); }
	if (dir != -1) Y = Inter(Y, Union(Yd, Yg));
	if (dir != 1) X = Inter(X, Union(Xd, Xg));
}
//----------------------------------------------------------------------
void Csqrt(interval& Y, interval& X, int dir)
{
	// Y=sqrt(X)           =>  dir=1;
	// X=sqrt-1(Y)         =>  dir=-1; 
	if (dir != -1) Y = Inter(Y, Sqrt(X));
	if (dir != 1) Csqr(X, Y, -dir);
}
//----------------------------------------------------------------------
void Cexp(interval& Y, interval& X, int dir)
{
	// Y=Exp(X)           =>  dir=1;
	// X=Exp-1(Y)=Log(Y)  =>  dir=-1; 
	if (dir != -1) { Y = Inter(Y, Exp(X)); }
	if (dir != 1) { X = Inter(X, Log(Y)); }
}
//----------------------------------------------------------------------
void Clog(interval& Y, interval& X, int dir)
{
	// Y=log(X)           =>  dir=1;
	// X=log-1(Y)=Exp(Y)  =>  dir=-1; 
	if (dir != -1) { Y = Inter(Y, Log(X)); }
	if (dir != 1) { X = Inter(X, Exp(Y)); }
}
//----------------------------------------------------------------------
void Cpow(interval& Y, interval& X, int n, int dir)
{
	// Y=Power(X,n)        =>   dir=1;
	// X=Power(Y,1/n)      =>   dir=-1; 
	if (n > 0)
	{
		interval Yd, Yg;
		interval Xd = Inter(X, interval(0, oo)), Xg = Inter(X, interval(-oo, 0));
		if (!Xd.isEmpty) 
		{
			Xd = Inter(Xd, PowRoot(Y, 1, n));
			Yd = Inter(Y, Pow(Xd, n));
		}
		if (!Xg.isEmpty) 
		{
			if (n % 2 == 0) Xg = Inter(Xg, -PowRoot(Y, 1, n));
			else Xg = Inter(Xg, PowRoot(Y, 1, n));
			Yg = Inter(Y, Pow(Xg, n));
		}
		if (dir != -1) Y = Inter(Y, Union(Yd, Yg));
		if (dir != 1) X = Inter(X, Union(Xd, Xg));
	}
	else 
	{
		interval Z = Pow(X, abs(n)), I(1);
		if (dir != -1) { Cpow(Z, X, abs(n), 1); Cdiv(Y, I, Z, 1); }
		if (dir != 1) { Cdiv(Y, I, Z, -1); Cpow(Z, X, abs(n), -1); }
	}
}
//----------------------------------------------------------------------
void Ccos(interval& Y, interval& X, int dir)
{
	// Y=Cos(X)    =>  dir=1;
	// X=Cos-1(Y)  =>  dir=-1; 
	if (dir != -1) { Y = Inter(Y, Cos(X)); }
	if (dir != 1) { interval X0 = X + PI_2; Csin(Y, X0, -1); X = Inter(X, X0 - PI_2); }
}
//-----------------------------------------------------------------------
void Csin(interval& Y, interval& X, int dir)
{
	// Y=sin(X)    =>  dir=1;
	// X=sin-1(Y)  =>  dir=-1;
	if (dir != -1) { Y = Inter(Y, Sin(X)); }
	if (dir != 1)
	{
		Y = Inter(Y, interval(-1, 1));
		if (Y.isEmpty) X = interval();
		else {
			interval Y1 = Inter(Inflate(Y, 1e-5), interval(-1, 1));
			double amoins = -oo, aplus = -oo,
				hmoins = asin(Y1.inf), hplus = asin(Y1.sup);
			int k;

			if (In(sin(X.inf), Y1)) { amoins = X.inf; }
			else  { // calcul de a-
				k = (int)(ceil((X.inf - hmoins) / DEUX_PI));
				amoins = hmoins + (k*DEUX_PI);
				k = (int)(ceil((X.inf + hmoins - PI) / DEUX_PI));
				amoins = min(amoins, PI - hmoins + (k*DEUX_PI));
				k = (int)(ceil((X.inf - hplus) / DEUX_PI));
				amoins = min(amoins, hplus + (k*DEUX_PI));
				k = (int)(ceil((X.inf + hplus - PI) / DEUX_PI));
				amoins = min(amoins, PI - hplus + (k*DEUX_PI));
			}
			if (amoins > X.sup) { X = interval(); }
			else {
				if (In(sin(X.sup), Y1)) { aplus = X.sup; }
				else {//calcul de a+
					k = (int)(floor((X.sup - hmoins) / DEUX_PI));
					aplus = hmoins + (k*DEUX_PI);
					k = (int)(floor((X.sup + hmoins - PI) / DEUX_PI));
					aplus = max(aplus, PI - hmoins + (k*DEUX_PI));
					k = (int)(floor((X.sup - hplus) / DEUX_PI));
					aplus = max(aplus, hplus + (k*DEUX_PI));
					k = (int)(floor((X.sup + hplus - PI) / DEUX_PI));
					aplus = max(aplus, PI - hplus + (k*DEUX_PI));
				}
			}
			X = Inter(X, interval(amoins, aplus));
		}
	}
}
//----------------------------------------------------------------------
void Ctan(interval& Y, interval& X, int dir)
{
	// Y=Tan(X)            =>  dir=1;
	// X=Tan-1(Y)=Atan(Y)  =>  dir=-1; 
	if (dir != -1) { Y = Inter(Y, Tan(X)); }
	if (dir != 1)
	{
		double amoins = -oo, aplus = -oo, hmoins = atan(Y.inf), hplus = atan(Y.sup);

		if (In(tan(X.inf), Y)) amoins = X.inf;
		else amoins = hmoins + (ceil((X.inf - hmoins) / PI)*PI);

		if (amoins > X.sup) X = interval();
		else 
		{
			if (In(tan(X.sup), Y)) aplus = X.sup;
			else aplus = hplus + (floor((X.sup - hplus) / PI)*PI);
		}
		X = Inter(X, interval(amoins, aplus));
	}
} 
//----------------------------------------------------------------------
void Catan(interval& Y, interval& X, int dir)
{
	// Y=atan(X)                  =>  dir=1;
	// X=atan-1(Y)=tan(Y)         =>  dir=-1; 
	Ctan(X, Y, -dir);
}
//----------------------------------------------------------------------
void Csinc(interval& Y, interval& X, int dir)
{
	// Y=sinc(X)=sin(X)/X      =>  dir=1;
	// X=sinc-1(X)             =>  dir=-1; 
	interval Z(-oo, oo);
	if (dir != -1) 
	{
		interval U(-oo, oo), V(-oo, oo), W(-oo, oo), G(-oo, oo);
		interval M(-oo, oo), N(-oo, oo), H(-oo, oo), I(-oo, oo);
		interval C(-oo, oo), un(1);
		/*Cabs(M,X,1); Cdiv(N,un,M,1);
		if (!N.isEmpty) Y=Inter(Y,interval(-oo,N.sup));
		Cabs(C,X,1); Csub(H,C,1); Cdiv(I,un,H,1);
		if (!I.isEmpty) Y=Inter(Y,interval(I.inf,oo));*/
		Csqr(G, X, 1); U = Inter(U, G / 6);
		Csub(V, U, 1); Cadd(W, V, un, 1);
		Y = Inter(Y, interval(-0.3, 1));
		if (!W.isEmpty) Y = Inter(Y, interval(W.inf, oo));
		Csin(Z, X, 1); Cdiv(Y, Z, X, 1);
	}
	if (dir != 1) { Cdiv(Y, Z, X, -1); Csin(Z, X, -1); }
}
//----------------------------------------------------------------------
void Carg(interval& A, interval& X, interval& Y, int dir)
{
	// A=arg(X,Y)                  =>  dir=1;
	// Y=X*tan(A);X=Y/tan(A);      =>  dir=-1; 
	if (dir != -1) 
	{
		//A.Intersect(Arg(X,Y));
		interval D = Sqrt(Sqr(X)+Sqr(Y)), y_x = Y/X, x_d = X/D, y_d = Y/D;
		Ctan(y_x, A, -1); Ccos(x_d, A, -1); Csin(y_d, A, -1);
	}
	if (dir != 1) 
	{
		interval y_tana = Y/Tan(A), xtana = X*Tan(A);
		X.Intersect(y_tana); Y.Intersect(xtana);
		interval D = Sqrt(Sqr(X)+Sqr(Y)), dcosa = D*Cos(A), dsina = D*Sin(A);
		X.Intersect(dcosa); Y.Intersect(dsina);
	}
}
//----------------------------------------------------------------------
void CAngle1(interval& X2, interval& Y2, interval& Theta, interval& X1, interval& Y1)
{
	interval SinTheta = Sin(Theta);
	interval CosTheta = Cos(Theta);
	interval CosThetaX1 = CosTheta*X1;
	interval SinThetaX1 = SinTheta*X1;
	interval CosThetaY1 = CosTheta*Y1;
	interval SinThetaY1 = SinTheta*Y1;
	X2 = Inter(CosThetaX1 - SinThetaY1, X2);
	Y2 = Inter(SinThetaX1 + CosThetaY1, Y2);
	Cadd(Y2, SinThetaX1, CosThetaY1, -1);
	Csub(X2, CosThetaX1, SinThetaY1, -1);
	Cmul(SinThetaY1, SinTheta, Y1, -1);
	Cmul(CosThetaY1, CosTheta, Y1, -1);
	Cmul(SinThetaX1, SinTheta, X1, -1);
	Cmul(CosThetaX1, CosTheta, X1, -1);
	Ccos(CosTheta, Theta, -1);
	Csin(SinTheta, Theta, -1);
}
//----------------------------------------------------------------------
int CAngle(interval& X2, interval& Y2, interval& Theta, interval& X1, interval& Y1, bool StrongAngle)
{
	CAngle1(X2, Y2, Theta, X1, Y1);
	if (StrongAngle)
	{
		interval Theta2 = -Theta;
		CAngle1(X1, Y1, Theta2, X2, Y2);
		Theta = -Theta2;
		interval SqrX1, SqrY1, SqrX2, SqrY2, N2;
		N2 = Inter(Sqr(X1) + Sqr(Y1), Sqr(X2) + Sqr(Y2));
		Cnorm(N2, X1, Y1);
		Cnorm(N2, X2, Y2);
	}
	if (X2.isEmpty || Y2.isEmpty || Theta.isEmpty || X1.isEmpty || Y1.isEmpty) return -1;
	return 1;
}
//----------------------------------------------------------------------
void Cnorm(interval& N, interval& X, interval& Y)
{
	interval SqrX, SqrY, SqrN;
	SqrX = Sqr(X);
	SqrY = Sqr(Y);
	SqrN = Sqr(N);
	Cadd(SqrN, SqrX, SqrY, 1);
	Cadd(SqrN, SqrX, SqrY, -1);
	Csqr(SqrY, Y, -1);
	Csqr(SqrX, X, -1);
	Csqr(SqrN, N, -1);
}
//----------------------------------------------------------------------
void Cnorm(interval& N, interval& X, interval& Y, interval& Z, int dir)
{
	interval SqrX = Sqr(X), SqrY = Sqr(Y), SqrZ = Sqr(Z), SqrN = Sqr(N);
	if (dir != -1)
	{
		SqrN = Inter(SqrN, SqrX+SqrY+SqrZ);
		Csqr(SqrN, N, -1);
	}
	if (dir != 1)
	{
		SqrZ = Inter(SqrZ, SqrN-SqrX-SqrY);
		SqrY = Inter(SqrY, SqrN-SqrX-SqrZ);
		SqrX = Inter(SqrX, SqrN-SqrZ-SqrY);
		Csqr(SqrZ, Z, -1);
		Csqr(SqrY, Y, -1);
		Csqr(SqrX, X, -1);
	}
}
//----------------------------------------------------------------------
void Cdist(interval& R, interval& X1, interval& Y1, interval& X2, interval& Y2)
{
	interval deltaX = X2 - X1; interval deltaY = Y2 - Y1;
	Cnorm(R, deltaX, deltaY);
	Csub(deltaX, X2, X1, -1); Csub(deltaY, Y2, Y1, -1);
}
//----------------------------------------------------------------------
void Cscal(interval& s, interval& ux, interval& uy, interval& vx, interval& vy)
{
	interval z1 = ux*vx;
	interval z2 = uy*vy;
	Cadd(s, z1, z2);
	Cmul(z2, uy, vy, -1);
	Cmul(z1, ux, vx, -1);
}
//----------------------------------------------------------------------
void Cscal(interval& s, double& ux, double& uy, interval& vx, interval& vy)
{
	interval z1 = ux*vx;
	interval z2 = uy*vy;
	Cadd(s, z1, z2);
	Cmul(z2, uy, vy, -1);
	Cmul(z1, ux, vx, -1);
}
//----------------------------------------------------------------------
void Cdet(interval& det, interval& ux, interval& uy, interval& vx, interval& vy, int dir)
{
	interval z1 = ux*vy;
	interval z2 = vx*uy;
	Csub(det, z1, z2, 1);
	if (dir == 1) return;
	Csub(det, z1, z2, -1);
	Cmul(z2, vx, uy, -1);
	Cmul(z1, ux, vy, -1);
}
//----------------------------------------------------------------------
void Cdet(interval& det, double& ux, double& uy, interval& vx, interval& vy, int dir)
{
	interval z1 = ux*vy;
	interval z2 = uy*vx;
	Csub(det, z1, z2, -1);
	if (dir == 1) return;
	Csub(det, z1, z2, -1);
	Cmul(z2, uy, vx, -1);
	Cmul(z1, ux, vy, -1);
}
//----------------------------------------------------------------------
void Cdet(interval& det, interval& ux, interval& uy, double& vx, double& vy, int dir)
{
	interval z1 = vy*ux;
	interval z2 = vx*uy;
	Csub(det, z1, z2, -1);
	if (dir == 1) return;
	Csub(det, z1, z2, -1);
	Cmul(z2, vx, uy, -1);
	Cmul(z1, vy, ux, -1);
}
//----------------------------------------------------------------------
void Cstep(interval& Y, interval& X)
{   
	vector<double> ax,ay,bx,by;
	ax.push_back(-1E10); ay.push_back(0); bx.push_back(0);    by.push_back(0);
	ax.push_back(0);     ay.push_back(0); bx.push_back(0);    by.push_back(1);
	ax.push_back(0);     ay.push_back(1); bx.push_back(1E10); by.push_back(1);
	CPointInSegments(X,Y,ax,ay,bx,by);
}
//----------------------------------------------------------------------
void Cstep(interval& Y, interval& X, int dir, double a)
{
	/* Y=step(X)=1 if X.inf>a
	=0 if X.sup<a
	=[0,1] if X.inf=<a<=X.sup   =>  dir=1

	X=step-1(Y)=Empty if Inter(Y,[0,1])=Empty
	=Empty if X.sup<a & 0 \not in Y
	=Empty if X.inf>a & 1 \not in Y
	=X if X.sup<a & 0 \in Y
	=X if X.inf>a & 1 \in Y
	=X if X.inf=<a<=X.sup   =>  dir=-1 */
	Y = Inter(Y, interval(0, 1));
	if (dir != -1)
	{
		if (!X.isEmpty)
		{
			if (X.inf > a) Y = Inter(Y, interval(1));
			if (X.sup < a) Y = Inter(Y, interval(0));
		}
		else Y = interval();
	}
	if (dir != 1)
	{
		if (!Y.isEmpty)
		{
			interval Yh(-oo, oo), Yb(-oo, oo), Xg(-oo, oo), Xd(-oo, oo);
			Yh = Inter(Y, interval(1)); Yb = Inter(Y, interval(0));
			Xg = Inter(X, interval(-oo, a)); Xd = Inter(X, interval(a, oo));
			if (Yh.isEmpty) Xd = Inter(Xd, interval(-oo, a));
			if (Yb.isEmpty) Xg = Inter(Xg, interval(a, oo));
			X = Inter(X, Union(Xg, Xd));
		}
		else X = interval();
	}
}
//----------------------------------------------------------------------
void Cramp(interval& Y, interval& X, int dir, double a)
{
	// Y=ramp(X)=max(0,X);
	// X=ramp-1(Y)
	Y = Inter(Y, interval(0, oo));
	if (dir != -1) { interval Zero(a); Cmax(Y, X, Zero, 1); }
	if (dir != 1) 
	{
		interval Xd = Inter(X, interval(a, oo));
		interval Xg = Inter(X, interval(-oo, a));
		if (Y.isEmpty) Xg = Xd = interval();
		else 
		{
			Xd = Inter(Xd, Y);
			if (Y.inf > 0) Xg = interval();
		}
		X = Union(Xd, Xg);
	}
}
//----------------------------------------------------------------------
void Cheaviside(interval& Y, interval& X, int dir, double a)
{
	/* Y=heavyside(X)=1 if X.inf>=a
	=0 if X.sup<a
	X=heavyside-1(Y)=Empty if Inter(Y,[0,1])=Empty
	=Empty if X.sup<a & 0 \not in Y
	=Empty if X.inf>a & 1 \not in Y
	=X if X.sup<a & 0 \in Y
	=X if X.inf>a & 1 \in Y
	=X if X.inf=<a<=X.sup   =>  dir=-1 */
	interval Z = Inter(interval(0), Y);
	interval W = Inter(interval(1), Y);
	if (dir != -1)
	{
		//interval U = Union(Z, W);
		if (X.isEmpty) Y = X;
		else 
		{
			if (Z.isEmpty && W.isEmpty) Y = Z;
			else 
			{
				Y = Inter(Y, Union(Z, W));
				if (X.inf >= a)  Y = Inter(Y, interval(1));
				else if (X.sup < a) Y = Inter(Y, interval(0));
			}
		}
	}
	if (dir != 1)
	{
		if (Z.isEmpty && W.isEmpty) X = Y = Z;
		else {
			if ((!Z.isEmpty) && (!W.isEmpty)) X = Inter(X, interval(-oo, oo));
			if (Z.isEmpty) X = Inter(X, interval(0, oo));
			if (W.isEmpty) X = Inter(X, interval(-oo, 0));
		}
	}
}
//----------------------------------------------------------------------
void Crect(interval& Z, interval& X, interval& Y, int dir)
{
	interval A1(-oo, oo), A2(-oo, oo), A3(-oo, oo), A4(-oo, oo), deux(2);
	X = Inter(X, interval(0, oo));
	interval aire = Y / deux;
	if (dir != -1)
	{
		A1 = Inter(A1, X + aire); Cstep(A2, A1, 1);
		A3 = Inter(A3, X - aire); Cstep(A4, A3, 1);
		Z = Inter(Z, A2 - A4);
	}
	if (dir != 1)
	{
		A2 = Inter(A2, interval(0, 1)); A4 = Inter(A4, interval(0, 1));
		A2 = Inter(A2, Z + A4); A4 = Inter(A4, A2 - Z);
		Cstep(A2, A1, -1); Cadd(A1, X, aire, -1); Y = Inter(Y, deux*aire);
		aire = Inter(aire, Y / deux); Cstep(A4, A3, -1); Csub(A3, X, aire, -1);
		Y = Inter(Y, deux*aire);
	}
}
//----------------------------------------------------------------------
void Crect(interval& Y, interval& X, int dir)
{
	interval V(1); Crect(Y, X, V, dir);
}
//----------------------------------------------------------------------
void Ctriangle(interval& Y, interval& X, int dir)
{
	interval Xg = Inter(X, interval(-oo, -0.5));
	interval Xd = Inter(X, interval(0.5, oo));
	interval Xm = Inter(X, interval(-0.5, 0.5));
	interval un(1), deux(2);
	Y = Inter(Y, interval(0, 1));
	if (dir != -1)
	{
		interval Yg(0), Yd(0), Ym(0);
		if (Xg.isEmpty) Yg = interval();
		if (Xd.isEmpty) Yd = interval();
		if (Xm.isEmpty) Ym = interval();
		else Ym = Inter(Y, (-deux*Abs(X)) + un);
		Y = Inter(Y, Union(Yd, Union(Yg, Ym)));
	}
	if (dir != 1)
	{
		if (Y.isEmpty) X = interval();
		else {
			interval A(-oo, oo);
			if (Y.inf > 0) Xg = Xd = interval();
			A = Inter(A, (un - Y) / deux);
			Cabs(A, Xm, -1);
			X = Inter(X, Union(Union(Xg, Xm), Xd));
		}
	}
}
//----------------------------------------------------------------------
void CDistanceDirLine(interval& dist, interval& mx, interval& my, interval& theta, double& ax, double& ay, double& bx, double& by)
{     
	// la distance dist entre le point m=(mx,my) a la droite [a,b] suivant le vecteur u
	if ((dist.isEmpty)||(mx.isEmpty)||(my.isEmpty)||(theta.isEmpty))
	{  
		dist = interval(); mx = interval(); my = interval(); theta = interval(); return;
	}
	interval ma_x=ax-mx; interval ma_y=ay-my;
	interval mb_x=bx-mx; interval mb_y=by-my;
	double ab_x=bx-ax; double ab_y=by-ay;
	interval ux=Cos(theta); interval uy=Sin(theta);
	interval z3=Det(ma_x,ma_y,ab_x,ab_y);
	interval z4=Det(ux,uy,ab_x,ab_y);
	dist=Inter(dist,z3/z4);

	Cdiv(dist,z3,z4, -1);
	Cdet(z4,ux,uy,ab_x,ab_y,-1);
	Cdet(z3,ma_x,ma_y,ab_x,ab_y,-1);
	Csin(uy,theta,-1); Ccos(ux,theta,-1);
	Csub(mb_y,by,my,-1); Csub(mb_x,bx,mx,-1);
	Csub(ma_y,ay,my,-1); Csub(ma_x,ax,mx,-1);
}
//----------------------------------------------------------------------
int CDistanceDirSegment(interval& dist, interval& mx, interval& my, interval& theta, double ax, double ay, double bx, double by, int dir)
{
	// la distance dist entre le point m=(mx,my) au segment [a,b] suivant le vecteur u
	if ((dist.isEmpty)||(mx.isEmpty)||(my.isEmpty)||(theta.isEmpty))
	{
		dist = interval(); mx = interval(); my = interval(); theta = interval(); return -1;
		//dist = interval(); mx = interval(); my = interval(); theta = interval(); return;
	}
	interval ma_x = ax - mx; interval ma_y = ay - my;
	interval mb_x = bx - mx; interval mb_y = by - my;
	double ab_x = bx - ax; double ab_y = by - ay;
	interval ux = Cos(theta); interval uy = Sin(theta);
	interval z1 = Det(ma_x, ma_y, ux, uy);
	interval z2 = Det(ux, uy, mb_x, mb_y);
	interval z3 = Det(ma_x, ma_y, ab_x, ab_y);
	interval z4 = Det(ux, uy, ab_x, ab_y);
	interval z5 = Min(z1, z2, z3);
	interval d1 = z3 / z4;
	interval infty = interval(oo, oo);
	Cchi(dist, z5, infty, d1);
	//REDONDANT
	Cchi(dist, z3, infty, dist);
	Cchi(dist, z4, infty, dist);
	if (dir == 1) return 1;
	//if (dir == 1) return;
	Cdiv(d1, z3, z4, -1);
	Cmin(z5, z1, z2, z3, -1);
	Cdet(z4, ux, uy, ab_x, ab_y, -1);
	Cdet(z3, ma_x, ma_y, ab_x, ab_y, -1);
	Cdet(z2, ux, uy, mb_x, mb_y, -1);
	Cdet(z1, ma_x, ma_y, ux, uy, -1);
	Csin(uy, theta, -1); Ccos(ux, theta, -1);
	Csub(mb_y, by, my, -1); Csub(mb_x, bx, mx, -1);
	Csub(ma_y, ay, my, -1); Csub(ma_x, ax, mx, -1);
	return 1;
}
//----------------------------------------------------------------------
void CDistanceDirSegments(interval& distmin, interval& mx, interval& my, interval& theta, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{      
	// Distance directionnelle relativement a un polygone
	vector<interval> dist(ax.size());
	for (unsigned int j = 0; j < ax.size(); j++) dist[j] = interval(0, oo);
	for (unsigned int j = 0; j < ax.size(); j++)
		CDistanceDirSegment(dist[j], mx, my, theta, ax[j], ay[j], bx[j], by[j], 1);
	Cmin(distmin, dist, -1);
	for (int j = (int)(ax.size() - 1); j >= 0; j--)
		CDistanceDirSegment(dist[j], mx, my, theta, ax[j], ay[j], bx[j], by[j], -1);
}
//----------------------------------------------------------------------
/*// Version normalement plus efficace
void CDistanceDirSegments(interval& distmin, interval& mx, interval& my, interval& theta, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{      
// Distance directionnelle relativement a un polygone
box P(3);   P[1]=mx;  P[2]=my;   P[3]=theta;

vector<interval> dist(ax.size());
for (int j=0;j<ax.size();j++) dist[j]=interval(0,oo);
vector <box> L(ax.size());
for (int j=0;j<ax.size();j++)
{  
L[j]=P;
CDistanceDirSegment(dist[j],P[1],P[2],P[3],ax[j],ay[j],bx[j],by[j],1);
}
//Cmin(distmin,dist,-1);
for (int j=0;j<ax.size();j++)
{  
vector <interval> Y(ax.size());
for (int j1=0;j1<ax.size();j1++)
if (j1==j) Y[j1]=distmin;
else Y[j1]=interval(distmin.inf,oo);
CDistanceDirSegment(dist[j],L[j][1],L[j][2],L[j][3],ax[j],ay[j],bx[j],by[j],1);
}
C_q_in(P,1,L);
mx=P[1];  my=P[2];   theta=P[3];
}*/
//----------------------------------------------------------------------
void CPointInLine(interval& mx, interval& my, double& ax, double& ay, double& bx, double& by)
{    
	// contracte relativement a la contrainte : "m appartient a la droite (a,b)"
	interval ma_x=ax-mx;
	interval ma_y=ay-my;
	double ab_x=bx-ax;
	double ab_y=by-ay;
	interval z1=interval(0,0);
	Cdet(z1,ab_x,ab_y,ma_x,ma_y,-1);
	Csub(ma_y,ay,my,-1);
	Csub(ma_x,ax,mx,-1);
	if ((mx.isEmpty)||(my.isEmpty)) { mx = interval(); my = interval(); }
}
//----------------------------------------------------------------------
void CPointInSegment(interval& mx, interval& my, double ax, double ay, double bx, double by)
{      
	// contracte relativement a la contrainte : "m appartient au segment [a,b]"
	mx = Inter(mx, Union(interval(ax, ax), interval(bx, bx)));
	my = Inter(my, Union(interval(ay, ay), interval(by, by)));
	interval ma_x = ax - mx;
	interval ma_y = ay - my;
	double ab_x = bx - ax;
	double ab_y = by - ay;
	interval z1 = interval(0, 0);
	Cdet(z1, ab_x, ab_y, ma_x, ma_y, -1);
	Csub(ma_y, ay, my, -1);
	Csub(ma_x, ax, mx, -1);
	if ((mx.isEmpty)||(my.isEmpty)) { mx = interval(); my = interval(); }
}
//----------------------------------------------------------------------
void CPointInSegments(interval& mx, interval& my, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{      
	// contracte relativement a la contrainte : "m appartient au polygone dont les segments sont les [ai,bi]"
	if ((ax.size() == 0)||(ay.size() == 0)||(bx.size() == 0)||(by.size() == 0)) return;
	vector<interval> Mx(ax.size());
	vector<interval> My(ax.size());
	for (unsigned int j = 0; j < ax.size(); j++)
	{
		interval mx0 = mx;
		interval my0 = my;
		CPointInSegment(mx0, my0, ax[j], ay[j], bx[j], by[j]);
		Mx[j] = mx0;
		My[j] = my0;
	}
	mx = Inter(mx, Union(Mx));
	my = Inter(my, Union(My));
}
//----------------------------------------------------------------------
void CPointInCircle(interval& mx, interval& my, double cx, double cy, double r)
{      
	// contracte relativement a la contrainte : "m appartient au cercle de centre c et de rayon r"
	interval d_x = mx - cx;
	interval d_y = my - cy;
	interval d2_x = Sqr(d_x);
	interval d2_y = Sqr(d_y);
	interval r2 = r*r;
	Cadd(r2, d2_x, d2_y, -1);
	Csqr(d2_x, d_x, -1);
	Csqr(d2_y, d_y, -1);
	mx = Inter(mx, d_x + cx);
	my = Inter(my, d_y + cy);
	if ((mx.isEmpty)||(my.isEmpty)) { mx = interval(); my = interval(); }
}
//----------------------------------------------------------------------
void CPointInCircles(interval& mx, interval& my, vector<double> cx, vector<double> cy, vector<double> r, bool truth)
{      
	// contracte relativement a la contrainte : "m appartient a un des cercles de centre ci et de rayon ri"
	if ((cx.size() == 0)||(cy.size() == 0)||(r.size() == 0)) return;
	vector<interval> Mx(cx.size());
	vector<interval> My(cx.size());
	for (unsigned int j = 0; j < cx.size(); j++)
	{
		interval mx0 = mx;
		interval my0 = my;
		CPointInCircle(mx0, my0, cx[j], cy[j], r[j]);
		Mx[j] = mx0;
		My[j] = my0;
	}
	if (truth)
	{
		mx = Inter(mx, Union(Mx));
		my = Inter(my, Union(My));
	} 
	else 
	{
		mx = Inter(mx, Inter(Mx));
		my = Inter(my, Inter(My));
	}
}
//----------------------------------------------------------------------
void CPointInSegmentsOrCircles(interval& mx, interval& my, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, vector<double> cx, vector<double> cy, vector<double> r)
{      
	// contracte relativement a la contrainte : "m appartient soit au polygone soit a un des cercles de centre ci et de rayon ri"
	interval mx1 = mx;
	interval my1 = my;
	CPointInSegments(mx1, my1, ax, ay, bx, by);
	interval mx2 = mx;
	interval my2 = my;
	CPointInCircles(mx2, my2, cx, cy, r);
	if ((ax.size() > 0)&&(cx.size() > 0))            //cercles et polygone
	{
		mx = Inter(mx, Union(mx1, mx2)); 
		my = Inter(my, Union(my1, my2));
	}
	if ((ax.size() > 0)&&(cx.size() == 0))           // pas de cercles
	{
		mx = mx1; my = my1;
	}
	if (ax.size() == 0)          // pas de segments
	{
		mx = mx2; my = my2;
	}
}
//----------------------------------------------------------------------
void CPointOutsideSegment(interval& mx, interval& my, double& ax, double& ay, double& bx, double& by, bool outer)
{
	// contracte relativement a la contrainte : "m appartient au segment [a,b]"
	//interval Ix = Union(interval(ax,ax),interval(bx,bx));
	//interval Iy = Union(interval(ay,ay),interval(by,by));
	interval I1x(-oo,oo), I1y(-oo,oo);

	//double dx = bx - ax; double dy = by - ay;
	//if(dx > 0 && dy > 0) { I1x = interval(ax,oo); I1y = interval(-oo, by); };
	//if(dx < 0 && dy > 0) { I1x = interval(bx,oo); I1y = interval(ay,oo); };
	//if(dx < 0 && dy < 0) { I1x = interval(-oo,ax); I1y = interval(by,oo); };
	//if(dx > 0 && dy < 0) { I1x = interval(-oo,bx); I1y = interval(-oo,ay); };

	mx = Inter(mx,I1x);
	my = Inter(my,I1y);
	//       mx=Inter(mx, Union(Ix,I1x));
	//       my=Inter(my, Union(Iy,I1y));
	//mx=Union(Inter(mx,Ix),Inter(mx,I1x));
	//my=Union(Inter(my,Iy),Inter(my,I1y));

	//       mx=Inter(mx,Union(interval(ax,ax),interval(bx,bx)));
	//       my=Inter(my,Union(interval(ay,ay),interval(by,by)));

	interval ma_x=ax-mx;
	interval ma_y=ay-my;
	interval mb_x=bx-mx;
	interval mb_y=by-my;
	interval z1;
	z1  = (outer == false) ? interval(0,oo) : interval(-oo,0);
	Cdet(z1,mb_x,mb_y,ma_x,ma_y,-1);
	Csub(ma_y,ay,my,-1);
	Csub(ma_x,ax,mx,-1);
	if ((mx.isEmpty)||(my.isEmpty)) { mx = interval(); my = interval(); }
}
//----------------------------------------------------------------------
void CPointOutsideSegments(interval& mx, interval& my, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, bool outer)
{ 
	if ((ax.size() == 0)||(ay.size() == 0)||(bx.size() == 0)||(by.size() == 0)) return;
	// contracte relativement a la contrainte : "m appartient au polygone dont les segments sont les [ai,bi]"
	vector<interval> Mx(ax.size());
	vector<interval> My(ax.size());
	for (unsigned int j = 0; j < ax.size(); j++)
	{  
		interval mx0=mx;
		interval my0=my;
		CPointOutsideSegment(mx0,my0,ax[j],ay[j],bx[j],by[j],false);
		Mx[j]=mx0;
		My[j]=my0;
	}
	if (outer == true)
	{
		mx=Inter(mx,Union(Mx));
		my=Inter(my,Union(My));
	}
	else
	{
		mx=Inter(mx,Inter(Mx));
		my=Inter(my,Inter(My));
	}
}
//----------------------------------------------------------------------
void CPoseInSegment(interval& mx, interval& my, interval& phi, double& ax, double& ay, double& bx, double& by)
{     
	// contracte relativement a "la pose (m,phi) appartient au segment [a,b]"
	CPointInSegment(mx,my,ax,ay,bx,by);
	double ab_x=bx-ax;   //(bx-ax)*cos(phi)+(by-ay)*sin(phi)=0
	double ab_y=by-ay;
	interval cphi=Cos(phi);
	interval sphi=Sin(phi);
	interval scal=interval(-0.0,0.0);   CScal(scal,ab_x,ab_y,cphi,sphi); // phi orthogonal to the line
	interval det(0,oo);
	Cdet(det,cphi,sphi,ab_x,ab_y);   // select the right direction
	Ccos(cphi,phi,-1);
	Csin(sphi,phi,-1);
	if ((mx.isEmpty)||(my.isEmpty)||(phi.isEmpty)) { mx = interval(); my = interval(); phi = interval(); }
}
//----------------------------------------------------------------------
void CPoseInSegments(interval& mx, interval& my, interval& phi, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{
	if ((ax.size() == 0)||(ay.size() == 0)||(bx.size() == 0)||(by.size() == 0)) return;
	vector<interval> Mx(ax.size());
	vector<interval> My(ax.size());
	vector<interval> Mphi(ax.size());
	for (unsigned int j=0;j<ax.size();j++)
	{ 
		interval mx0=mx;
		interval my0=my;
		interval phi0=phi;
		CPoseInSegment(mx0,my0,phi0,ax[j],ay[j],bx[j],by[j]);
		Mx[j]=mx0;
		My[j]=my0;
		Mphi[j]=phi0;
	}
	mx=Inter(mx,Union(Mx));
	my=Inter(my,Union(My));
	phi=Inter(phi,Union(Mphi));
}
//----------------------------------------------------------------------
void CPoseInCircle(interval& mx, interval& my, interval& phi, double& cx, double& cy, double& r)
{ 
	CPointInCircle(mx,my,cx,cy,r);
	interval mc_x=cx-mx;
	interval mc_y=cy-my;
	interval cphi=Cos(phi);
	interval sphi=Sin(phi);
	interval scal=interval(0,oo);   CScal(scal,mc_x,mc_y,cphi,sphi); //scal(cphi,sphi,cx_mx,cy-my)>0
	interval det(0,0);
	Cdet(det,cphi,sphi,mc_x,mc_y);   //det(cphi,sphi,cx_mx,cy-my)=0
	Ccos(cphi,phi,-1);
	Csin(sphi,phi,-1);
	if ((mx.isEmpty)||(my.isEmpty)||(phi.isEmpty)) { mx = interval(); my = interval(); phi = interval(); }
}
//----------------------------------------------------------------------
void CPoseInCircles(interval& mx, interval& my, interval& phi,vector<double> cx, vector<double> cy, vector<double> r)
{ 
	if ((cx.size() == 0)||(cy.size() == 0)||(r.size() == 0)) return;
	vector<interval> Mx(cx.size());
	vector<interval> My(cx.size());
	vector<interval> Mphi(cx.size());
	for (unsigned int j=0;j<cx.size();j++)
	{  
		interval mx0=mx;
		interval my0=my;
		interval phi0=phi;
		CPoseInCircle(mx0,my0,phi0,cx[j],cy[j],r[j]);
		Mx[j]=mx0;
		My[j]=my0;
		Mphi[j]=phi0;
	}
	mx=Inter(mx,Union(Mx));
	my=Inter(my,Union(My));
	phi=Inter(phi,Union(Mphi));
}
//----------------------------------------------------------------------
void CPoseInSegmentsOrCircles(interval& mx, interval& my, interval& malpha, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by,
							  vector<double> cx, vector<double> cy, vector<double> r)
{      
	// "la pose (m,alpha) appartient soit au polygone soit a un des cercles de centre ci et de rayon ri"
	vector<interval> Mx(2);
	vector<interval> My(2);
	vector<interval> Malpha(2);
	interval mx0=mx;
	interval my0=my;
	interval malpha0=malpha;
	CPoseInSegments(mx0,my0,malpha0,ax,ay,bx,by);
	Mx[0]=mx0;   My[0]=my0;  Malpha[0]=malpha0;
	mx0=mx;  my0=my;  malpha0=malpha;
	CPoseInCircles(mx0,my0,malpha0,cx,cy,r);
	Mx[1]=mx0;
	My[1]=my0;
	Malpha[1]=malpha0;
	mx=Inter(mx,Union(Mx));
	my=Inter(my,Union(My));
	malpha=Inter(malpha,Union(Malpha));
	// Bug if no segments or no circles, see CPointInSegmentsOrCircles()...?
}
//----------------------------------------------------------------------
void CPoseTrans(interval& qx, interval& qy, interval& dist, interval& px, interval& py, interval& alpha)
{ 
	interval ux=Cos(alpha);
	interval uy=Sin(alpha);
	interval dx=ux*dist;
	interval dy=uy*dist;
	interval pxdx = px+dx;
	interval pydy = py+dy;
	qx=Inter(qx,pxdx);
	qy=Inter(qy,pydy);
	Cadd(qy,py,dy,-1);
	Cadd(qx,px,dx,-1);
	Cmul(dy,uy,dist,-1);
	Cmul(dx,ux,dist,-1); // Plante ici car dx est negatif et d > 0 et ux > 0?
	Csin(uy,alpha,-1);
	Ccos(ux,alpha,-1);
	if ((qx.isEmpty)||(qy.isEmpty)||(alpha.isEmpty)) { qx = interval(); qy = interval(); alpha = interval(); }
}
//----------------------------------------------------------------------
void CPoseTransRot(interval& qx, interval& qy, interval& beta, interval& d, interval& psi, interval& px,interval& py,interval& alpha)
{  
	CPoseTrans(qx,qy,d, px, py,alpha);
	Cadd(beta,psi,alpha);
}
//----------------------------------------------------------------------
void CPoseRotTrans(interval& qx, interval& qy, interval& beta, interval& phi,interval& d, interval& px,interval& py,interval& alpha)
{   
	Cadd(beta,phi,alpha);
	CPoseTrans(qx,qy,d, px, py,beta);
	Cadd(beta,phi,alpha);
}
//----------------------------------------------------------------------
void CPoseRotTransRot(interval& qx, interval& qy, interval& beta, interval& phi, interval& d, interval& psi, interval& px,interval& py,interval& alpha)
{   
	interval delta=alpha+phi;
	CPoseTransRot(qx,qy,beta,d,psi,px,py,delta);
	Cadd(delta,alpha,phi,-1);
}
//----------------------------------------------------------------------
void CPoseTransInWallsOrCircles(interval& px, interval& py, interval& alpha, interval& d,
								vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by,
								vector<double> cx, vector<double> cy, vector<double> r)
{    
	interval qx=interval(-oo,oo);
	interval qy=interval(-oo,oo);
	CPoseTrans(qx,qy,d,px, py,alpha);
	CPoseInSegmentsOrCircles(qx, qy, alpha, ax, ay, bx, by,cx,cy,r);
	CPoseTrans(qx,qy,d,px,py,alpha);
}
//----------------------------------------------------------------------
void CPoseTransRotInWallsOrCircles(interval& px, interval& py, interval& alpha, interval& d, interval& psi, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, vector<double> cx, vector<double> cy, vector<double> r)
{   
	// All feet should be supported by one of the circles or one of the segments
	interval qx=interval(-oo,oo);
	interval qy=interval(-oo,oo);
	interval beta=interval(-oo,oo);
	CPoseTransRot(qx,qy,beta,d, psi,px, py,alpha);
	CPoseInSegmentsOrCircles(qx, qy, beta, ax, ay, bx, by,cx,cy,r);
	CPoseTransRot(qx,qy,beta,d, psi,px, py,alpha);
}
//----------------------------------------------------------------------
void CPoseRotTransRotInWallsOrCircles(interval& px, interval& py, interval& alpha, interval& phi, interval& d, interval& psi, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, vector<double> cx, vector<double> cy, vector<double> r)
{  
	interval qx=interval(-oo,oo);
	interval qy=interval(-oo,oo);
	interval beta=interval(-oo,oo);
	CPoseRotTransRot(qx,qy,beta,phi,d,psi,px,py,alpha);
	CPoseInSegmentsOrCircles(qx, qy, beta, ax, ay, bx, by,cx,cy,r);
	CPoseRotTransRot(qx,qy,beta,phi,d,psi,px,py,alpha);
}
//----------------------------------------------------------------------
void CPoseRotTransPointInWallsOrCircles(interval& px, interval& py, interval& alpha, interval& phi, interval& d, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, vector<double> cx, vector<double> cy, vector<double> r)
{    
	interval qx=interval(-oo,oo);
	interval qy=interval(-oo,oo);
	interval beta=interval(-oo,oo);
	CPoseRotTrans(qx,qy,beta,phi,d,px,py,alpha);
	CPointInSegmentsOrCircles(qx, qy, ax, ay, bx, by,cx,cy,r);
	CPoseRotTrans(qx,qy,beta,phi,d,px,py,alpha);
}
//----------------------------------------------------------------------
void CPoseTransPointInLine(interval& px, interval& py, interval& alpha, interval& d, double ax, double ay, double bx, double by)
{    
	interval qx=interval(-oo,oo);
	interval qy=interval(-oo,oo);
	CPoseTrans(qx,qy,d, px, py,alpha);
	CDistanceDirLine(d,px, py, alpha,ax, ay,bx,by);
	CPointInLine(qx,qy,ax,ay,bx,by);
	CPoseTrans(qx,qy,d, px, py,alpha);
}
//----------------------------------------------------------------------
void CPoseTransPointInCircles(interval& px, interval& py, interval& alpha, interval& d, vector<double>& cx, vector<double>& cy, vector<double>& r, bool truth)
{   
	interval qx=interval(-oo,oo);
	interval qy=interval(-oo,oo);
	CPoseTrans(qx,qy,d, px, py,alpha);
	CPointInCircles(qx,qy,cx,cy,r,truth);
	CPoseTrans(qx,qy,d, px, py,alpha);
}
//----------------------------------------------------------------------
void CPoseTransPointInWall(interval& px, interval& py, interval& alpha, interval& d0, double ax, double ay, double bx, double by, bool truth)
{    
	//interval qx=interval(-oo,oo);
	//interval qy=interval(-oo,oo);
	interval d=interval(-oo,oo);
	interval px1(px);
	interval py1(py);
	interval alpha1(alpha);
	interval px2(px);
	interval py2(py);
	interval alpha2(alpha);
	//qDebug()<<"px="<<px<<"py="<<py<<"alpha="<<alpha<<"d="<<d;

	CPoseTransPointInLine(px1,py1,alpha1,d,ax,ay,bx,by);
	if (truth) d=Inter(d,d0); else Cnotin(d,d0);
	CPoseTransPointInLine(px1,py1,alpha1,d,ax,ay,bx,by);

	CPoseTowardSegment(px2, py2, alpha2, ax, ay,bx,by,truth);

	if (truth) {px=Inter(px1,px2);py=Inter(py1,py2);alpha=Inter(alpha1,alpha2);}  // Intersection des deux contracteurs
	else {px=Union(px1,px2);py=Union(py1,py2);alpha=Union(alpha1,alpha2);}        // Union des deux contracteurs

}
//----------------------------------------------------------------------
void CPoseTransPointInWalls(interval& px, interval& py, interval& alpha, interval& d0, vector<double>& ax, vector<double>& ay, vector<double>& bx, vector<double>& by, bool truth)
{
	vector<interval> Lpx;
	vector<interval> Lpy;

	for(unsigned int i = 0; i < ax.size(); i++){
		interval px1(px); interval py1(py); interval d(d0); interval theta(alpha);
		double eps = 0.25;
		double dx = (bx[i] - ax[i]);
		double dy = (by[i] - ay[i]);
		double dxn = (dx/hypot(dx,dy))*eps;
		double dyn = (dy/hypot(dx,dy))*eps;

		//CPoseTransPointInWall(px1,py1,theta,d,ax[i],ay[i],bx[i],by[i],truth);
		CPoseTransPointInWall(px1,py1,theta,d,ax[i]-dxn,ay[i]-dyn,bx[i]+dxn,by[i]+dyn,truth);
		Lpx.push_back(px1);
		Lpy.push_back(py1);
	}
	if (truth)
	{
		px = Inter(px, Union(Lpx));
		py = Inter(py, Union(Lpy));
	} 
	else 
	{
		px = Inter(px, Inter(Lpx));
		py = Inter(py, Inter(Lpy));
	}
}
//----------------------------------------------------------------------
void CPoseTransPointInWallsOrCircles(interval& px, interval& py, interval& alpha, interval& d0, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by,
									 vector<double> cx, vector<double> cy, vector<double> r, bool truth)
{
	vector<interval> Px(2), Py(2), Theta(2);
	interval px0, py0, theta0;
	if (truth)
	{
		interval d(d0);
		//px0 = px; py0 = py, theta0 = alpha; d=d0;
		//CPoseTransPointInWalls(px0,py0,theta0,d,ax,ay,bx,by,truth);
		//Px[0] = px0; Py[0]=py0; Theta[0] = theta0;

		px0 = px; py0 = py, theta0 = alpha; d = d0;
		CPoseTransPointInCircles(px,py,theta0,d,cx,cy,r, truth);
		Px[0] = px0; Py[0] = py0; Theta[0] = theta0;

		//px = Inter(px, Union(Px));
		//py = Inter(py, Union(Py));
		//alpha = Inter(alpha, Union(Theta));
	} 
	else 
	{
		//px0 = px; py0 = py, theta0 = alpha ;
		//CPoseTransPointInWalls(px0,py0,theta0,d0,ax,ay,bx,by,truth);
		//Px[0] = px0; Py[0] = py0; Theta[0] = theta0;

		px0 = px; py0 = py, theta0 = alpha;
		interval d(0,oo);
		CPoseTransPointInCircles(px0,py0,theta0,d, cx,cy,r,truth);
		Cnotin(d,d0);
		CPoseTransPointInCircles(px,py,theta0,d, cx,cy,r,truth);
		//Px[1] = px0; Py[1] = py0; Theta[1] = theta0;
		//px = Inter(px, Inter(Px));
		//py = Inter(py, Inter(Py));
		//alpha = Inter(alpha, Inter(Theta));
	}
}
//----------------------------------------------------------------------
void CPoseTowardSegment(interval& mx, interval& my, interval& theta, double& ax, double& ay, double& bx,double& by, bool truth)
{      
	// La pose m=(mx,my,theta) pointe sur le segment [a,b]
	if ((mx.isEmpty)||(my.isEmpty)||(theta.isEmpty))
	{ 
		mx = interval(); my = interval(); theta = interval(); return;
	}
	interval ma_x=ax-mx;       interval ma_y=ay-my;
	interval mb_x=bx-mx;       interval mb_y=by-my;
	double ab_x=bx-ax;         double ab_y=by-ay;
	interval ux=Cos(theta);    interval uy=Sin(theta);
	interval z1=Det(ma_x,ma_y,ux,uy);
	interval z2=Det(ux,uy,mb_x,mb_y);
	interval z3=Det(ma_x,ma_y,ab_x,ab_y);
	interval z4=Det(ux,uy,ab_x,ab_y);

	interval z5=interval(0,oo);      // si tous les zi >0 alors on satisfait la contrainte
	if (!truth) z5=interval(-oo,0);
	Cmin(z5,z1,z2,z3,z4,-1);

	//z1=Inter(z1,interval(0,oo));
	//z2=Inter(z2,interval(0,oo));
	//z3=Inter(z3,interval(0,oo));
	//z4=Inter(z4,interval(0,oo));

	Cdet(z4,ux,uy,ab_x,ab_y,-1);
	Cdet(z3,ma_x,ma_y,ab_x,ab_y,-1);
	Cdet(z2,ux,uy,mb_x,mb_y,-1);
	Cdet(z1,ma_x,ma_y,ux,uy,-1);
	Csin(uy,theta,-1);           Ccos(ux,theta,-1);
	Csub(mb_y,by,my,-1);       Csub(mb_x,bx,mx,-1);
	Csub(ma_y,ay,my,-1);       Csub(ma_x,ax,mx,-1);
}
//----------------------------------------------------------------------
void Cnocross(interval& px, interval& py, interval& mx, interval& my, double& ax, double& ay, double& bx, double& by)
{      
	interval ma_x=ax-mx;
	interval ma_y=ay-my;
	interval mb_x=bx-mx;
	interval mb_y=by-my;
	interval ap_x=px-ax;
	interval ap_y=py-ay;
	interval pb_x=bx-px;
	interval pb_y=by-py;
	interval pm_x=mx-px;
	interval pm_y=my-py;
	double ab_x=bx-ax;
	double ab_y=by-ay;
	interval z1=interval(-oo,oo);
	Cdet(z1,ab_x,ab_y,ma_x,ma_y,1);
	interval z2=interval(-oo,oo);
	Cdet(z2,ab_x,ab_y,ap_x,ap_y,1);
	interval z3=z1*z2;
	interval z4=interval(-oo,oo);
	Cdet(z4,pm_x,pm_y,ap_x,ap_y,1);
	interval z5=interval(-oo,oo);
	Cdet(z5,pm_x,pm_y,pb_x,pb_y,1);
	interval z6=z4*z5;
	interval z7(0,oo);
	Cmax(z7,z3,z6,-1);
	Cmul(z6,z4,z5,-1);
	Cdet(z5,pm_x,pm_y,pb_x,pb_y,-1);
	Cdet(z4,pm_x,pm_y,ap_x,ap_y,-1);
	Cmul(z3,z1,z2,-1);
	Cdet(z2,ab_x,ab_y,ap_x,ap_y,-1);
	Cdet(z1,ab_x,ab_y,ma_x,ma_y,-1);
	Csub(pm_y,my,py,-1);
	Csub(pm_x,mx,px,-1);
	Csub(pb_y,by,py,-1);
	Csub(pb_x,bx,px,-1);
	Csub(ap_y,py,ay,-1);
	Csub(ap_x,px,ax,-1);
	Csub(mb_y,by,my,-1);
	Csub(mb_x,bx,mx,-1);
	Csub(ma_y,ay,my,-1);
	Csub(ma_x,ax,mx,-1);
	if ((mx.isEmpty)||(my.isEmpty)||(px.isEmpty)||(py.isEmpty))
	{  
		mx = interval(); my = interval(); px = interval(); py = interval(); 
	}
}
//----------------------------------------------------------------------
void CLegCrossNoSegment(interval& dist, interval& px, interval& py, interval& theta, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{      
	// Aucun segment ne doit être croise
	if ((ax.size() == 0)||(ay.size() == 0)||(bx.size() == 0)||(by.size() == 0)) return;
	interval ux=Cos(theta);
	interval uy=Sin(theta);
	interval dx=ux*dist;
	interval dy=uy*dist;
	interval leg_x=px+dx;
	interval leg_y=py+dy;
	for (unsigned int j = 0; j < ax.size(); j++)
	{
		Cnocross(px,py,leg_x,leg_y, ax[j], ay[j], bx[j], by[j]);
	}
	Cadd(leg_y,py,dy,-1);
	Cadd(leg_x,px,dx,-1);
	Cmul(dy,uy,dist,-1);
	Cmul(dx,ux,dist,-1);
	Csin(uy,theta,-1);
	Ccos(ux,theta,-1);
}
//----------------------------------------------------------------------
void CLegOnWalls(interval& dist, interval& px, interval& py, interval& theta, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by)
{      
	// Toutes les pattes doivent être sur le mur
	if ((ax.size() == 0)||(ay.size() == 0)||(bx.size() == 0)||(by.size() == 0)) return;
	interval ux = Cos(theta);
	interval uy = Sin(theta);
	interval dx = ux*dist;
	interval dy = uy*dist;
	interval leg_x = px + dx;
	interval leg_y = py + dy;
	CPointInSegments(leg_x, leg_y, ax, ay, bx, by);
	Cadd(leg_y, py, dy, -1);
	Cadd(leg_x, px, dx, -1);
	Cmul(dy, uy, dist, -1);
	Cmul(dx, ux, dist, -1);
	Csin(uy, theta, -1);
	Ccos(ux, theta, -1);
}
//----------------------------------------------------------------------
void CLegOnWallsOrCircles(interval& dist, interval& px, interval& py, interval& theta, vector<double> ax, vector<double> ay, vector<double> bx, vector<double> by, vector<double> cx, vector<double> cy, vector<double> r)
{      
	// Toutes les pattes doivent être sur le mur ou sur un des cercles
	interval ux = Cos(theta);
	interval uy = Sin(theta);
	interval dx = ux*dist;
	interval dy = uy*dist;
	interval leg_x = px + dx;
	interval leg_y = py + dy;
	CPointInSegmentsOrCircles(leg_x, leg_y, ax, ay, bx, by, cx, cy, r);
	Cadd(leg_y, py, dy, -1);
	Cadd(leg_x, px, dx, -1);
	Cmul(dy, uy, dist, -1);
	Cmul(dx, ux, dist, -1);
	Csin(uy, theta, -1);
	Ccos(ux, theta, -1);
}
//-----------------------------------------------------------------------
//             PROCEDURES DE CONTRACTION ELEMENTAIRE
//-----------------------------------------------------------------------
/* Ceci est le code choisit pour designes chacune des fontions elementaires
ABS    <=>   'B'   ?
COS    <=>   'C'
SIN    <=>   'S'
TAN    <=>   'T'
ATAN   <=>   'A'     ?
EXP    <=>   'E'
LOG    <=>   'L'      ?
SQR    <=>   '2'
SQRT   <=>   'R'
MAX    <=>   'M'       ?
MIN    <=>   'm'       ?
POWER  <=>   'P'
-Unaire <=>  '_'
Integer <=>  'z'        ?
Boolean <=>  'b'        ?
door
rect
step
sinC
*/
//-----------------------------------------------------------------------
void Contract0(char code, interval& Y)
{
	if (code == 'b')    { Cboolean(Y);    return; }
	if (code == 'z')    { Cinteger(Y);   return; }
}
//-----------------------------------------------------------------------
void Contract0(char code, interval& Y, interval& X, int dir)
	/* Procedure de propagation de contraintes binaires sur les
	intervalles, cette procedure est chargee de faire appel a une
	des fonctions binaires.
	exemple:  Y=sin(X)     =>   Contract0('S',Y,X,1)
	X=sin-1(Y)   =>   Contract0('S',Y,X,-1) */

{
	if (code == 'S')   { Csin(Y, X, dir);     return; }
	if (code == 'C')   { Ccos(Y, X, dir);     return; }
	if (code == 'E')   { Cexp(Y, X, dir);     return; }
	if (code == 'T')   { Ctan(Y, X, dir);     return; }
	if (code == '2')   { Csqr(Y, X, dir);     return; }
	if (code == 'R')   { Csqrt(Y, X, dir);    return; }
	if (code == 'L')   { Clog(Y, X, dir);     return; }
	if (code == 'A')   { Catan(Y, X, dir);    return; }
	if (code == '_')   { Csub(Y, X, dir);   return; }
	if (code == '=')   { Cequal(Y, X, dir);    return; }
	if (code == 'B')   { Cabs(Y, X, dir);     return; }
	if (code == 'D')   { Crect(Y, X, dir);    return; }
	if (code == 's')   { Csinc(Y, X, dir);     return; }
	//if (code=='B')   {Cpow(Y,X,2,dir); return;}
}
//----------------------------------------------------------------------
/*void Contract0 (char code, interval& Y, interval& X, int n, int dir)
Procedure de propagation de contraintes concernant la puissance entiere
sur les intervalles : Y=X^n, cette procedure est chargee de faire appel a
l'operateur elementaire Power.
exemple:  Y=Power(X,n)     =>   Contract0('P',Y,X,n,1)
X=Power-1(Y,n)   =>   Contract0('P',Y,X,n,-1)
{ if (code=='P') Cpow(Y,X,n,dir); } */
//----------------------------------------------------------------------
void Contract0(char code, interval& Z, interval& Y, interval& X, int dir)
	/* Procedure de propagation de contraintes ternaires sur les
	intervalles, cette procedure est chargee de faire appel aux
	fonctions ternaires.
	exemple:  Z=Y+X         =>  Contract0('+',Z,Y,X,1)
	Y=Z-X; X=Z-Y  =>  Contract0('+',Z,Y,X,-1) */
{
	if (code == '+')    { Cadd(Z, Y, X, dir);    return; }
	if (code == '-')    { Csub(Z, Y, X, dir);   return; }
	if (code == '*')    { Cmul(Z, Y, X, dir);    return; }
	if (code == '/')    { Cdiv(Z, Y, X, dir);     return; }
	if (code == 'm')    { Cmin(Z, Y, X, dir);     return; }
	if (code == 'M')    { Cmax(Z, Y, X, dir);     return; }
}
//----------------------------------------------------------------------
void ShowContraction(interval& Xcd, interval& Xcg, interval& X, interval& Xc)
{
	if (Xc.isEmpty) { Xcd = Xcg = Xc; return; }
	Xcd = interval(X.inf, Xc.inf);
	Xcg = interval(Xc.sup, X.sup);
}
//----------------------------------------------------------------------
void IntButterfly(interval& Y, interval Yo, interval dY, interval& X, interval Xo, int dir)
{
	UNREFERENCED_PARAMETER(dir);

	interval z1(-oo, +oo);
	interval z2(-oo, +oo);

	Csub(z1, X, Xo, 1);
	Cmul(z2, dY, z1, 1);
	Cadd(Y, Yo, z2, 1);

	Cadd(Y, Yo, z2, -1);
	Cmul(z2, dY, z1, -1);
	Csub(z1, X, Xo, -1);
}
//----------------------------------------------------------------------
void Inter1(interval& r0, interval& r1, const interval &a, const interval &b, const interval &c)
{   
	//interval r(-oo);
	if (a.isEmpty) { r0 = a; r1 = Inter(b, c); return; };
	if (b.isEmpty) { r0 = b; r1 = Inter(a, c); return; };
	if (c.isEmpty) { r0 = c; r1 = Inter(a, b); return; };
	double zmin0 = NAN, zmax0 = NAN, zmin1 = NAN, zmax1 = NAN;
	//Min1(zmin0,zmin1,Inf(a),Inf(b),Inf(c));
	//Max1(zmax0,zmax1,Sup(a),Sup(b),Sup(c));
	Max1(zmin0, zmin1, a.inf, b.inf, c.inf);
	Min1(zmax0, zmax1, a.sup, b.sup, c.sup);
	if (zmin0 > zmax0) r0 = interval(); else r0 = interval(zmin1, zmax1);
	if (zmin1 > zmax1) r1 = interval(); else r1 = interval(zmin1, zmax1);
	return;
}
//----------------------------------------------------------------------
void Sucre(interval& P, const interval& S)
{
	if (Disjoint(P, S)||Subset(S, P)||Subset(P, S)) return;
	if (In(S.inf, P)) { P.sup = S.inf; return; }
	if (In(S.sup, P)) { P.inf = S.sup; return; }
}
//----------------------------------------------------------------------
//void Cnotin(interval& X, const interval& Y, int dir)
//{  
//	UNREFERENCED_PARAMETER(dir);
//	if (Y.isEmpty) return;
//	interval X1 = Inter(X, interval(-oo,Y.inf));
//    interval X2 = Inter(X, interval(Y.sup,+oo));
//    X = Inter(X, Union(X1,X2));
//}
//----------------------------------------------------------------------
void Cnotin(interval& X, interval& Y, int dir)
{ 
	UNREFERENCED_PARAMETER(dir);
	if ((X.isEmpty)||(Y.isEmpty)) return;
	else
	{
		if (X.inf <= Y.inf)
		{
			if ((X.sup <= Y.sup)&&(Y.inf <= X.sup)) X = interval(X.inf, Y.inf);
			else return;
		}
		else
		{
			if (X.sup <= Y.sup) X = interval();
			else
			{
				if (X.inf <= Y.sup) X = interval(Y.sup, X.sup);
				else return;
			}
		}
	}
}
//----------------------------------------------------------------------
void C_q_in(interval& x, int q, vector<interval>& y)
{
	vector<borne> V;
	if (x.isEmpty) return; // Added.
	if (q <= 0) return; // Added.
	for (unsigned int i = 0; i < y.size(); i++)
	{
		if (!y[i].isEmpty) // Added.
		{
			V.push_back(borne(y[i].inf, 1));
			V.push_back(borne(y[i].sup, -1));
		}
	}
	if (V.size() == 0) // Added.
	{
		x = interval();
		return;
	}
	V.push_back(borne(x.inf, 2)); // Moved, removed -0.00001 and replaced 0 by 2.
	V.push_back(borne(x.sup, -2)); // Moved, removed +0.00001 and replaced 0 by -2.
	sort(V.begin(), V.end());
	int sum = 0;
	int imin = (int)(V.size());
	int imax = -1;
	for (unsigned int i = 0; i < V.size(); i++)
	{
		if (abs(V[i].ouverture) != 2) sum = sum + V[i].ouverture; // Added condition to avoid x contribution to the sum.
		if ((sum >= q)&&(V[i].ouverture > 0)&&(V[i].val >= x.inf)&&(V[i].val <= x.sup)) // Replaced all & by &&, == 1 by > 0.
		{
			imin = min(imin, (int)i);
		}
		if ((sum >= q - 1)&&(V[i].ouverture < 0)&&(V[i].val >= x.inf)&&(V[i].val <= x.sup)) // Replaced all & by &&, == -1 by < 0.
		{
			imax = max(imax, (int)i);
		}
	}
	if ((imax == -1)||(imin == (int)(V.size()))) x = interval(); // Added condition to check also imin.
	else x = interval(V[imin].val, V[imax].val);
}
//----------------------------------------------------------------------
int q_in(double x, vector<interval>& y)
{
	int n = 0;
	for (unsigned int i = 0; i < y.size(); i++)
		if (In(x, y[i])) n++;
	return n;
}
//----------------------------------------------------------------------
/*
void C_q_in (interval& x, int q, vector<interval>& y)
{   
vector<double> V;
V.push_back(x.inf);
V.push_back(x.sup);
for (int i=0;i<y.size();i++)
{V.push_back(y[i].inf);
V.push_back(y[i].sup);
};
vector<double> E;
for (int i=0;i<V.size();i++)
{ 
if (In(V[i],x) & (q_in(V[i],y)>=q))
E.push_back(V[i]);
}
x=Enveloppe(E);
} 
*/
/*
void C_q_in(interval& x, int q, vector<interval>& y)
{
vector<borne> V;
V.push_back(borne(x.inf - 0.00001, 0));
V.push_back(borne(x.sup + 0.00001, 0));
for (unsigned int i = 0; i < y.size(); i++)
{
V.push_back(borne(y[i].inf, 1));
V.push_back(borne(y[i].sup, -1));
};
sort(V.begin(), V.end());
int sum = 0;
int imin = (int)(V.size());
int imax = -1;
for (unsigned int i = 0; i < V.size(); i++)
{
sum = sum + V[i].ouverture;
if ((sum >= q)&(V[i].ouverture == 1)&(V[i].val >= x.inf)&(V[i].val <= x.sup))
{
imin = min(imin, (int)i);
};
if ((sum >= q - 1)&(V[i].ouverture == -1)&(V[i].val >= x.inf)&(V[i].val <= x.sup))
{
imax = max(imax, (int)i);
};
}
if (imax == -1) x = interval();
else x = interval(V[imin].val, V[imax].val);
}*/
//----------------------------------------------------------------------
// Other
//----------------------------------------------------------------------
//modifs???
// RETURN LA DIFFERENCE ENTRE 2 intervalles x0, x1 et place le resultat dans c0 et c1.
// ON suppose que x0 inclus dans x1.
void diffI(interval &x0, interval &x1, interval &c0, interval &c1)
{
	//interval xt = Inter(Inter(x0, x1)-x1,interval(0,0));
	if (x1.isEmpty)
	{
		c0 = interval(); c1 = interval();
	} 
	else 
	{
		c0 = (x1.inf == x0.inf) ? interval() : interval(x1.inf,x0.inf);
		c1 = (x0.sup == x1.sup) ? interval() : interval(x0.sup,x1.sup);
		if (abs(c0.sup - c0.inf) < 1e-10) c0 = interval();
		if (abs(c1.sup - c1.inf) < 1e-10) c1 = interval();
	}
	if (c0.isEmpty)
	{
		c0 = c1;
		c1 = interval();
	}
}
//----------------------------------------------------------------------
iboolean TestDiskExists(const interval& X,const interval& Y,const interval& P1,const interval& P2,const interval& P3)
{
	// Test for the constraint : Exists p1 in P1, p2  in P2, (x-p1)^2+(y-p2)^2 in P3
	interval A1=X-P1.inf;
	interval B1=X-P1.sup;
	interval A2=Y-P2.inf;
	interval B2=Y-P2.sup;
	interval A12=Sqr(A1);
	interval A22=Sqr(A2);
	interval B12=Sqr(B1);
	interval B22=Sqr(B2);
	interval Z1=Step(B1*A1)*Min(B12,A12)+Step(B2*A2)*Min(B22,A22);
	interval Z2=Max(B12,A12)+Max(B22,A22);
	return (In(Min(Z2,P3.sup*P3.sup)-Max(Z1,P3.inf*P3.inf),interval(0,oo)));
}
//----------------------------------------------------------------------
iboolean TestDiskForall(const interval& X,const interval& Y,const interval& P1,const interval& P2,const interval& P3)
{   
	// Test for the constraint : For all p1 in P1, p2  in P2, (x1-p1)^2+(x2-p2)^2 in P3
	iboolean dedans1,dedans2;
	interval A1=X-P1.inf;
	interval B1=X-P1.sup;
	interval A2=Y-P2.inf;
	interval B2=Y-P2.sup;
	interval A12=Sqr(A1);
	interval A22=Sqr(A2);
	interval B12=Sqr(B1);
	interval B22=Sqr(B2);
	interval Z1=Step(B1*A1)*Min(B12,A12)+Step(B2*A2)*Min(B22,A22);
	interval Z2=Max(B12,A12)+Max(B22,A22);
	dedans1=In(Z2-P3.sup*P3.sup,interval(-oo,0));
	dedans2=In(Z1-P3.inf*P3.inf,interval(0,oo));
	return (dedans1&&dedans2);
}
//----------------------------------------------------------------------

#ifdef __GNUC__
// Restore the GCC warnings previously disabled.
#if (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 6)) || (__GNUC__ > 4))
#pragma GCC diagnostic pop
#else
//#if (__GNUC__ >= 9)
//#pragma GCC diagnostic warning "-Wdeprecated-copy"
//#endif // (__GNUC__ >= 9)
#endif // (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 6)) || (__GNUC__ > 4))
#endif // __GNUC__
