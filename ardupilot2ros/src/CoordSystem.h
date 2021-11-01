/***************************************************************************************************************:')

CoordSystem.h

Coordinate system handling.

Fabrice Le Bars

Created : 2009-02-12

***************************************************************************************************************:)*/

#ifndef COORDSYSTEM_H
#define COORDSYSTEM_H

#include "OSCore.h"

/*
Structure that defines a coordinate system.
*/
struct _COORDSYSTEM
{
	double xMin;
	double xMax;
	double yMin;
	double yMax;
};
typedef struct _COORDSYSTEM COORDSYSTEM;

#endif // !COORDSYSTEM_H
