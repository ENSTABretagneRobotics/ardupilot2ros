// Simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS and Jeremy NICOLA.

#include "iboolean.h"

#ifdef __GNUC__
// Disable some GCC warnings.
#if (__GNUC__ >= 9)
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#endif // (__GNUC__ >= 9)
#if (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 6)) || (__GNUC__ > 4))
#pragma GCC diagnostic push
#endif // (((__GNUC__ == 4) && (__GNUC_MINOR__ >= 6)) || (__GNUC__ > 4))
#endif // __GNUC__

iboolean::iboolean()
{
	value = iperhaps;
}

iboolean::iboolean(bool b)
{
	if (b) value = itrue; else value = ifalse;
}

iboolean::iboolean(IBOOLEAN tt)
{
	value = tt;
}

iboolean::iboolean(const iboolean& t)
{
	*this = t;
}

iboolean operator&&(iboolean x, iboolean y)
{
	return And(x,y);
}

iboolean operator||(iboolean x, iboolean y)
{
	return Or(x,y);
}

bool operator==(iboolean x, iboolean y)
{
	return (x.value == y.value);
}

bool operator!=(iboolean x, iboolean y)
{
	return (x.value != y.value);
}

iboolean operator!(iboolean x)
{
	if (x.value == itrue) return iboolean(ifalse);
	if (x.value == ifalse) return iboolean(itrue);
	return x;
}

iboolean Not(iboolean x)
{	
	if (x.value == iperhaps) return iperhaps;
	if (x.value == itrue) return ifalse;
	if (x.value == ifalse) return itrue;
	return iempty;
}

iboolean Inter(iboolean x, iboolean y)
{
	if (x.value == y.value) return x.value;
	if (x.value == iperhaps) return y.value;
	if (y.value == iperhaps) return x.value;
	return iboolean(iempty);
}

iboolean Union(iboolean x, iboolean y)
{	
	if (x.value == iempty) return y.value;
	if (y.value == iempty) return x.value;
	if (x.value == y.value) return x.value;
	return iboolean(iperhaps);
}

iboolean Xor(iboolean x, iboolean y)
{	
	if (x.value == iempty) return iempty;
	if (y.value == iempty) return iempty;
	if (x.value == iperhaps) return iperhaps;
	if (y.value == iperhaps) return iperhaps;
	if (x.value == y.value) return ifalse;
	return iboolean(itrue);
}

iboolean And(iboolean x, iboolean y)
{	
	if ((x.value == iempty)||(y.value == iempty)) return iboolean(iempty);
	if ((x.value == ifalse)||(y.value == ifalse)) return iboolean(ifalse);
	if ((x.value == iperhaps)||(y.value == iperhaps)) return iboolean(iperhaps);
	return iboolean(itrue);
}

// Project the constraint x>=y with respect to x
iboolean geq(iboolean x, iboolean y)      
{	
	iboolean r;
	r = iboolean(iperhaps);
	if (y.value == iempty) r = iboolean(iempty);
	if (y.value == itrue) r = iboolean(itrue);
	return Inter(x,r);
}

// Project the constraint x<=y with respect to x
iboolean leq(iboolean x, iboolean y)      
{	
	iboolean r;
	r = iboolean(iperhaps);
	if (y.value == iempty)  r = iboolean(iempty);
	if (y.value == ifalse) r = iboolean(ifalse);
	return Inter(x,r);
	/*                               a     &&    (implique b)
	1*0=\iempty                       1     &&      0
	1*1=1                            1     &&     [0,1]
	1*[0,1]=1                        1     &&     [0,1]
	0*0=0                            0     &&      0
	0*1=0                            0     &&      [0,1]
	0*[0,1]=0                        0     &&      [0,1]
	[0,1]*0=0                        [0,1] &&       0
	[0,1]*1=[0,1]                    [0,1] &&      [0,1]
	[0,1]*[0,1]=[0,1]                [0,1] &&     [0,1]
	*/
}

// return x and !y
iboolean Restrict(iboolean x, iboolean y) 
{ 
	return And(x,!y);
}

iboolean Or(iboolean x, iboolean y)
{	
	if ((x.value == iempty)||(y.value == iempty)) return iboolean(iempty);
	if ((x.value == itrue)||(y.value == itrue)) return iboolean(itrue);
	if ((x.value == iperhaps)||(y.value == iperhaps)) return iboolean(iperhaps);
	return iboolean(ifalse);
}

std::ostream& operator<<(std::ostream& os, const iboolean& a)
{
	if (a.value == itrue) os << "itrue";
	if (a.value == ifalse) os << "ifalse";
	if (a.value == iperhaps) os << "iperhaps";
	if (a.value == iempty) os << "iempty";
	return os;
}

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
