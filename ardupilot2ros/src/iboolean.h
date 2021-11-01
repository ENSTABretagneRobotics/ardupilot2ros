// Simple interval library from Luc JAULIN, with minor modifications from Fabrice LE BARS and Jeremy NICOLA.

#ifndef IBOOLEAN_H
#define IBOOLEAN_H

// IBOOLEAN is an interval Boolean. Used for a trivalued logic
// if x=[0,0]=ifalse  : certainly false
// if x=[0,1]=iperhaps: don't know
// if x=[1,1]=itrue  : certainly true
// Otherwise x=iempty

#ifdef _MSC_VER
// Enable additional features in math.h.
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif // !_USE_MATH_DEFINES
#endif // _MSC_VER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cfloat>

enum IBOOLEAN {itrue, ifalse, iperhaps, iempty};

class iboolean
{
public:
	IBOOLEAN value;
public:
	iboolean();
	iboolean(bool);
	iboolean(IBOOLEAN);
	iboolean(const iboolean&);
	friend iboolean operator&&(iboolean, iboolean);
	friend iboolean operator||(iboolean, iboolean);
	friend iboolean operator!(iboolean);
	friend bool operator==(iboolean, iboolean);
	friend bool operator!=(iboolean, iboolean);
	friend iboolean Not(iboolean x);
	friend iboolean Inter(iboolean, iboolean);
	friend iboolean Union(iboolean, iboolean);
	friend iboolean And(iboolean,iboolean);
	friend iboolean Or(iboolean,iboolean);
	friend iboolean leq(iboolean x, iboolean y);
	friend iboolean geq(iboolean x, iboolean y);
	friend iboolean Restrict(iboolean x, iboolean y); // return x and !y
	friend iboolean Xor(iboolean x, iboolean y);
	friend std::ostream& operator<<(std::ostream& os, const iboolean&);
};

#endif // !IBOOLEAN_H
